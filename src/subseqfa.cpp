/*
 * subseqfa.cpp
 *
 * Extract sub-sequences from records in an alignment and write a FASTA.
 *
 *  Created on: Nov 13, 2017
 *      Author: paudano
 */

#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <regex>

#include "seqtools/constants.h"
#include "seqtools/util/err.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/algorithm/string/split.hpp"

#include "htslib/sam.h"

// Namespaces
using namespace std;

namespace po = boost::program_options;

// Declarations
char *progName;

const char SEQI_TO_CHAR[] {
	'*', 'A', 'C', '*',
	'G', '*', '*', '*',
	'T', '*', '*', '*',
	'*', '*', '*', 'N',
};

// Main
int main(int argc, char *argv[]) {

	// Declarations
	string region;                  // Region
	vector<string> inFileNameList;  // List of input files
	string outFileName;             // Output file name
	bool verbose;                   // Verbose flag
	bool base0Region;               // Region is in base-0 (BED) coordinates

	// Region
	string chr;
	int pos;
	int end;
	vector<string> regionTok;

	// Output
	ofstream *outStream;  // Output file stream
	streambuf *outBuf;    // Output buffer (from outStream or cout)

	outStream = nullptr;

	// Alignment Objects
	samFile *inFile;      // Alignment input file
	bam_hdr_t *inHeader;  // Header for input file

	bam1_t *alignRecord = bam_init1();  // Alignment record

	const uint32_t *cigar;  // Array of CIGAR operations
	int nCigar;      // Number of CIGAR operations
	int cigarIndex;  // Position in the array of CIGAR operations

	uint8_t *seq;  // Sequence

	// Objects for subsetting
	int nextPos;    // Set to pos or end; CIGAR is traversed until nextPos is found
	int subPos;     // Sub-sequence start position
	int subEnd;     // Sub-sequence end
	int subIndex;   // Index for traversing from subPos to subEnd

	int32_t posTemplate;  // Template (reference) position
	int32_t posQuery;     // Query (read) position

	uint32_t cigarOp;       // CIGAR operation code
	uint32_t cigarLen;      // CIGAR operation length


	// Set program name for error reporting
	progName = argv[0];

	// Get program options
	po::options_description prog_opts("Extract regions from alignments and write to a FASTA file.");

	prog_opts.add_options()
			("help,h", "Print help")
			("region,r", po::value<string>(&region), "Region to extract (1-based, inclusive, chr:start-end)")
			("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print verbose information")
			("out,o", po::value<string>(&outFileName)->default_value(""), "Table of phased reads and HET-SNV counts")
			("base0,b", po::bool_switch(&base0Region)->default_value(false), "Region is in base-0 half-open coordinates (BED coordinates)")
			;

	po::options_description hidden_opts("Hidden options");
		hidden_opts.add_options()
				("infile", po::value<vector<string>>(&inFileNameList), "Input alignment file(s)")
				;

	po::positional_options_description p;
	p.add("infile", -1);

	// Merge options
	po::options_description cmdline_opts;
	cmdline_opts.add(prog_opts).add(hidden_opts);

	po::options_description visible("Allowed options");
	visible.add(prog_opts);

	// Parse options
	po::variables_map vm;
	try {
		po::store(po::command_line_parser(argc, argv).options(cmdline_opts).positional(p).run(), vm);

		po::notify(vm);

	} catch (const std::exception &ex) {
		err(ex);
		return ERR_USAGE;
	}

	// Print help
	if (vm.count("help")) {
		cout << progName << " [<options>] input1.sam/bam/cram [input2...]\n" << endl;
		cout << prog_opts << endl;
		return ERR_NONE;
	}

	// Parse region string
	boost::trim(region);
	region = regex_replace(region, regex("\\s+"), "-");

	boost::split(regionTok, region, boost::is_any_of(":-_"));

	if (regionTok.size() != 3) {
		err("Malformed region: Expected chr:pos-end (where delimiters may be :, -, _, or whitespace): \"%s\"", region.c_str());
	}

	chr = regionTok[0];
	pos = stoi(regionTok[1]);
	end = stoi(regionTok[2]);

	if (verbose)
		cout << "Region: " << chr << ":" << pos << "-" << end << endl;

	if (! base0Region)
		pos -= 1;

	// Open output
	boost::trim(outFileName);

	if (outFileName == "") {
		outBuf = cout.rdbuf();
	} else {
		outStream = new ofstream(outFileName);
		outBuf = outStream->rdbuf();
	}

	ostream out(outBuf);

	// Read files
	for (string inFileName : inFileNameList) {
		if (verbose)
			cout << "Reading " << inFileName << endl;

		inFile = hts_open(inFileName.c_str(), "r");

		inHeader = sam_hdr_read(inFile);

		// Read records
		while (sam_read1(inFile, inHeader, alignRecord) >= 0) {

			if (verbose)
				cout << "Record: " << (char *) alignRecord->data << endl;

			// Check chr, pos, and end (region must be fully-contained within the read)
			if (chr != inHeader->target_name[alignRecord->core.tid]) {
				if (verbose)
					cout << "\t* No region for target: " << inHeader->target_name[alignRecord->core.tid] << endl;

				continue;
			}

			if (alignRecord->core.pos > pos || bam_endpos(alignRecord) < end) {

				if (verbose)
					cout << "\t* Record does not cover query region" << endl;

				continue;
			}

			// Get cigar operations
			cigar = bam_get_cigar(alignRecord);
			nCigar = alignRecord->core.n_cigar;
			cigarIndex = 0;

			posTemplate = alignRecord->core.pos;

			// Traverse CIGAR and locate positions
			nextPos = pos;
			posQuery = 0;

			while (cigarIndex < nCigar && nextPos >= 0) {
				cigarOp = cigar[cigarIndex] & BAM_CIGAR_MASK;
				cigarLen = cigar[cigarIndex] >> BAM_CIGAR_SHIFT;

				switch(cigarOp) {

				// Matching bases
				case BAM_CMATCH:
				case BAM_CEQUAL:
				case BAM_CDIFF:

					if (posTemplate + cigarLen > nextPos) {

						if (nextPos == pos) {
							subPos = posQuery + (nextPos - posTemplate);
							nextPos = end;
							continue;

						} else {
							subEnd = posQuery + (nextPos - posTemplate);
							nextPos = -1;
							continue;
						}
					}

					posQuery += cigarLen;
					posTemplate += cigarLen;

					break;

				// Inserted bases (skip)
				case BAM_CINS:

					posQuery += cigarLen;

					break;

				// Deleted and skipped bases
				case BAM_CDEL:
				case BAM_CREF_SKIP:

					if (posTemplate + cigarLen > nextPos) {
						if (nextPos == pos) {
							subPos = posQuery;
							nextPos = end;
							continue;

						} else {
							subEnd = posQuery;
							nextPos = -1;
							continue;
						}
					}

					posTemplate += cigarLen;

					break;

				// Clipped bases
				case BAM_CSOFT_CLIP:

					posQuery += cigarLen;

					break;

				// Ignore other operations (hard clipping, padding)
				}  // switch CIGAR op

				// Next CIGAR op
				cigarIndex += 1;

			}  // Loop cigar ops


			//
			// Write FASTA record
			//

			if (subEnd <= subPos) {
				cout << "\t* No sequence found: " << subPos << "-" << subEnd << endl;
				continue;
			}

			// Prepare
			if (verbose)
				cout << "\t* Extracting: " << (char *) alignRecord->data << ":" << subPos << "-" << subEnd << endl;

			seq = bam_get_seq(alignRecord);

			// Write header
			out << ">" << alignRecord->data << ":" << subPos << "-" << subEnd;

			for (subIndex = subPos; subIndex < subEnd; ++subIndex) {

				if ((subIndex - subPos) % 80 == 0)
					out << "\n";

				out << SEQI_TO_CHAR[bam_seqi(seq, subIndex)];
			}

			// End record
			if ((subIndex - subPos) % 80 == 0)
				out << flush;
			else
				out << endl;


		}  // Loop records
	}  // Loop alignment input files

	return ERR_NONE;
}

