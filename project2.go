package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

////////////////////////////////////////////////////////////////////////////////
// Data structure & struct declaration
////////////////////////////////////////////////////////////////////////////////

type spacer struct {
	spacers map[string]location
}

type location struct {
	strand     bool
	PAM        string
	Startpoint int
	chromosome string
}

type param struct {
	GCmin    int
	GCmax    int
	hairpin  int
	upstream int
	contU    int
	addtlG   bool
}

////////////////////////////////////////////////////////////////////////////////
// open & read files
////////////////////////////////////////////////////////////////////////////////

// return filelist within directory
func readDir() []string {
	var filelist []string
	files, _ := filepath.Glob("*")
	filelist = append(files)
	return filelist
}

// read each file line by line
// from mark.go assignment
func readFile(filename string) []string {
	// open the file and make sure all went well
	in, err := os.Open(filename)
	if err != nil {
		fmt.Println("Error: couldn’t open the file")
		os.Exit(3)
	}
	// create the variable to hold the lines
	var lines = make([]string, 0)
	// for every line in the file
	scanner := bufio.NewScanner(in)
	for scanner.Scan() {
		// append it to the lines slice
		lines = append(lines, scanner.Text())
	}
	// check that all went ok
	if scanner.Err() != nil {
		fmt.Println("Sorry: there was some kind of error during the file reading")
		os.Exit(3)
	}
	// close the file and return the lines
	in.Close()
	return lines
}

// read and update the parameter struct
func readParam(file string) *param {
	var p param
	parameters := readFile(file)
	p.GCmin, _ = strconv.Atoi(parameters[2][7:])
	p.GCmax, _ = strconv.Atoi(parameters[3][7:])
	p.hairpin, _ = strconv.Atoi(parameters[4][9:])
	p.upstream, _ = strconv.Atoi(parameters[5][10:])
	p.contU, _ = strconv.Atoi(parameters[6][7:])
	p.addtlG, _ = strconv.ParseBool(parameters[7][8:])
	return &p
}

// create a single string of each chromosome
func readChromosome(file string) map[string]string {
	in, err := os.Open(file)
	if err != nil {
		fmt.Println("Error: couldn’t open the file")
		os.Exit(3)
	}
	var allChromosomes map[string]string
	allChromosomes = make(map[string]string)
	scanner := bufio.NewScanner(in)
	chromosome := ""
	key := ""
	for scanner.Scan() {
		// append it to the lines slice
		if strings.Contains(scanner.Text(), ">") {
			allChromosomes[key] = chromosome
			key = scanner.Text()[1:]
			chromosome = ""
			continue
		}
		chromosome += scanner.Text()
	}
	allChromosomes[key] = chromosome
	if scanner.Err() != nil {
		fmt.Println("Sorry: there was some kind of error during the file reading")
		os.Exit(3)
	}
	in.Close()
	delete(allChromosomes, "")
	return allChromosomes
}

// initialize the spacer map of struct
func initSpacers() *spacer {
	s := spacer{
		spacers: map[string]location{
			"": location{true, "", 0, ""},
		},
	}
	return &s
}

////////////////////////////////////////////////////////////////////////////////
// spacer requirements / tests
////////////////////////////////////////////////////////////////////////////////

// check all the requirements as set in parameter for spacer
func (p *param) checkSpacer(seq string) bool {
	// check if spacer contains "N"
	if strings.Contains(strings.ToUpper(seq), "N") {
		return false
	}
	if p.checkGCcontent(seq) {
		return false
	}
	if p.multipleU(seq) {
		return false
	}
	if p.RNAfold(seq) {
		return false
	}
	return true
}

// check if spacer satisfies min and max GC content rule
// return true if GCcontent is too low or too high
func (p *param) checkGCcontent(seq string) bool {
	var count int
	var percentage float64
	seq = strings.ToUpper(seq)
	for i := 0; i < len(seq); i++ {
		if string(seq[i]) == "G" || string(seq[i]) == "C" {
			count++
		}
	}
	percentage = float64(count) / float64(len(seq))
	if percentage < float64(p.GCmin)/100 {
		return true
	} else if percentage > float64(p.GCmax)/100 {
		return true
	}
	return false
}

// checks if there are multiple 'U's in the gRNA
// this translates to multiple 'A's in the genome
// returns true if there are runs of 'U's as defined by parameter
func (p *param) multipleU(seq string) bool {
	var Ucount int
	if p.contU == 0 { // user has chosen to ignore this parameter
		return false
	}
	for i := 0; i < len(seq); i++ {
		if string(seq[i]) == "A" {
			Ucount++
			if Ucount == p.contU {
				return true
			}
		} else {
			Ucount = 0
		}
	}
	return false
}

// returns true if there are additional G's after PAM when not allowed
func (p *param) trailingG(base string) bool {
	if p.addtlG {
		return false
	}
	if base == "G" {
		return true
	}
	return false
}

// check if 2 bases are complementary
func isComplement(s, t string) bool {
	s = strings.ToUpper(s)
	t = strings.ToUpper(t)
	if s == "A" && t == "T" {
		return true
	} else if s == "T" && t == "A" {
		return true
	} else if s == "C" && t == "G" {
		return true
	} else if s == "G" && t == "C" {
		return true
	}
	return false
}

// recursively check the longest extension of hairpin
func checkFold(s string, i, j int) int {
	if i == 0 || j == len(s)-1 {
		return 0
	}
	if isComplement(string(s[i]), string(s[j])) {
		return 1 + checkFold(s, i-1, j+1)
	}
	return 0
}

// return the lnogest hairpin
func hairpin(s string) int {
	n := len(s)
	var longest, length int
	for i := 0; i < n-4; i++ {
		for j := 4; j < n; j++ {
			if j-i < 5 {
				continue
			}
			length = checkFold(s, i, j)
			if length > longest {
				longest = length
			}
		}
	}
	return longest
}

// return true if seq contains hairpins
func (p *param) RNAfold(seq string) bool {
	longest := hairpin(seq)
	if longest >= p.hairpin {
		return true
	}
	return false
}

////////////////////////////////////////////////////////////////////////////////
// find & update the spacers
////////////////////////////////////////////////////////////////////////////////

// go through entire chromosomes and update the spacer map
func (s *spacer) findSpacer(allChromosomes map[string]string, p *param) {
	var doubles string // look at 2 sequences at a time - search for GG or CC
	var strand bool
	var PAM string
	var seq string
	// go through each chromosome
	for chromosome, sequence := range allChromosomes {
		// go through each position within sequence
		for index := 0; index < len(sequence)-1; index++ {
			doubles = sequence[index : index+2]
			// check the leading strand for PAM
			if doubles == "GG" && index > 20 {
				seq = sequence[index-21 : index-1]
				if p.checkSpacer(seq) {
					// check if entry exists
					if p.trailingG(string(sequence[index+2])) {
						continue
					}
					if s.checkKey(seq, p) { // key exists
						PAM = "---"
					} else { // key doesn't exist
						PAM = sequence[index-1 : index+2]
						strand = true
					}
					// update spacer
					s.update(seq, location{strand, PAM, index, chromosome})
				}
				// check the reverse complemenetary strand for PAM
			} else if doubles == "CC" && index < len(sequence)-22 {
				seq = oppositeStrand(sequence[index+3 : index+23])
				if p.checkSpacer(seq) {
					if p.trailingG(oppositeStrand(string(sequence[index-1]))) {
						continue
					}
					if s.checkKey(seq, p) { // key exists
						PAM = "---"
					} else { // key doesn't exist
						PAM = oppositeStrand(sequence[index : index+3])
						strand = false
					}
					// update spacer
					s.update(seq, location{strand, PAM, index, chromosome})
				}
			}
		}
	}
}

// return reverse complementary strand (in 3'-5') of input DNA sequence
func reverseComplement(seq string) string {
	var ComplementaryStrand string
	for base := 0; base < len(seq); base++ {
		switch {
		case strings.ToUpper(seq[base:base+1]) == "A":
			ComplementaryStrand += "T"
		case strings.ToUpper(seq[base:base+1]) == "C":
			ComplementaryStrand += "G"
		case strings.ToUpper(seq[base:base+1]) == "G":
			ComplementaryStrand += "C"
		case strings.ToUpper(seq[base:base+1]) == "T":
			ComplementaryStrand += "A"
		case strings.ToUpper(seq[base:base+1]) == "N":
			ComplementaryStrand += "N"
		}
	}
	return ComplementaryStrand
}

// return opposite strand (in 5'-3') of input DNA sequence
func oppositeStrand(seq string) string {
	var ComplementaryStrand string
	for base := len(seq) - 1; base >= 0; base-- {
		switch {
		case strings.ToUpper(seq[base:base+1]) == "A":
			ComplementaryStrand += "T"
		case strings.ToUpper(seq[base:base+1]) == "C":
			ComplementaryStrand += "G"
		case strings.ToUpper(seq[base:base+1]) == "G":
			ComplementaryStrand += "C"
		case strings.ToUpper(seq[base:base+1]) == "T":
			ComplementaryStrand += "A"
		case strings.ToUpper(seq[base:base+1]) == "N":
			ComplementaryStrand += "N"
		}
	}
	return ComplementaryStrand
}

// check if key exists in map, with allowed mismatches
func (s *spacer) checkKey(seq string, p *param) bool {
	for key := range s.spacers {
		if len(key) < 20 {
			continue
		} else {
			if key[len(key)-p.upstream-1:] == seq[len(seq)-p.upstream-1:] {
				return true
			}
		}
	}
	return false
}

// update the spacer map struct with new values
func (s *spacer) update(spacerseq string, data location) {
	s.spacers[spacerseq] = data
}

// remove duplicates
func (s *spacer) postProcess(p *param) {
	for key := range s.spacers {
		if s.spacers[key].PAM == "---" { // non-unique entries
			delete(s.spacers, key)
		} else if s.spacers[key].PAM == "" {
			delete(s.spacers, key)
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// output/print formatting
////////////////////////////////////////////////////////////////////////////////

// format parameter for output File
func (p *param) paramFormat() []string {
	var parameters []string
	parameters = append(parameters, "Parameters")
	gcmin := "Minimum GC content = " + strconv.Itoa(p.GCmin)
	parameters = append(parameters, gcmin)
	gcmax := "Maximum GC content = " + strconv.Itoa(p.GCmax)
	parameters = append(parameters, gcmax)
	hairpin := "Number of bases for secondary structure = " + strconv.Itoa(p.hairpin)
	parameters = append(parameters, hairpin)
	upstream := "Position at which mismatch is allowed = " + strconv.Itoa(p.upstream)
	parameters = append(parameters, upstream)
	contU := "maximum number of continuous U's allowed = " + strconv.Itoa(p.contU)
	parameters = append(parameters, contU)
	addtlG := "Additional G's after PAM site allowed? = " + strconv.FormatBool(p.addtlG)
	parameters = append(parameters, addtlG)
	parameters = append(parameters, "")
	parameters = append(parameters, "*****************************")
	parameters = append(parameters, "")
	return parameters
}

// outputs all spacers in user-friendly format
func (s *spacer) spacerFormat(key string) string {
	PAM := s.spacers[key].PAM
	chromosome := s.spacers[key].chromosome
	direction := strconv.FormatBool(s.spacers[key].strand)
	startlocus := strconv.Itoa(s.spacers[key].Startpoint)
	var endlocus string
	if s.spacers[key].strand {
		endlocus = strconv.Itoa(s.spacers[key].Startpoint + 20)
	} else {
		endlocus = strconv.Itoa(s.spacers[key].Startpoint - 20)
	}
	IDTReadySeq := IDTseq(key + PAM)
	return key + PAM + " | " + chromosome + " | " + direction + " | " + startlocus + " | " + endlocus + " | " + IDTReadySeq
}

// return IDT-order ready gRNA sequence
func IDTseq(target string) string {
	tracrRNA := "ACTTTTTCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCATAGCTCTTAAACATTTGTGTCCAAGAATGTTTCCCTATAGTGAGTCGTATTAACTTTTTCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCATAGCTCTTAAACATTTGTGTCCAAGAATGTTTCCCTATAGTGAGTCGTATTA"
	antisense := reverseComplement(target)
	return antisense + tracrRNA
}

////////////////////////////////////////////////////////////////////////////////
// print / write to file
////////////////////////////////////////////////////////////////////////////////

// find match between unique sites in the genome and the target sequence
func (s *spacer) findMatch(x *spacer, p *param) {
	outFile, err := os.Create("Target_Sequence.txt")
	if err != nil {
		fmt.Println("Error creating file")
	}
	// print out the settings/parameter used
	parameterFormatted := p.paramFormat()
	for _, line := range parameterFormatted {
		fmt.Fprintln(outFile, line)
	}
	for skey := range s.spacers {
		for xkey := range x.spacers {
			if skey == xkey {
				keyFormatted := s.spacerFormat(skey)
				fmt.Fprintln(outFile, keyFormatted)
			}
		}
	}
	outFile.Close()
}

// there is no target sequence
func (s *spacer) printAll(p *param) {
	outFile, err := os.Create("All_Target_Sites.txt")
	if err != nil {
		fmt.Println("Error creating file")
	}
	parameterFormatted := p.paramFormat()
	for _, line := range parameterFormatted {
		fmt.Fprintln(outFile, line)
	}
	for key := range s.spacers {
		keyFormatted := s.spacerFormat(key)
		fmt.Fprintln(outFile, keyFormatted)
	}
	outFile.Close()
}

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

func main() {
	fmt.Println("Running...")
	filelist := readDir()
	p := readParam("parameter.txt")
	s := initSpacers()      // for genome
	x := initSpacers()      // for target gene
	var t map[string]string // initialize t for target gene in case it is provided
	for _, file := range filelist {
		if strings.Contains(file, "chr") { // read chromosomes
			fmt.Println("Now reading : ", file)
			chromosomes := readChromosome(file)
			s.findSpacer(chromosomes, p)
		} else if file == "target.txt" { //read target
			fmt.Println("Now reading : target sequence")
			t = readChromosome("target.txt")
			x.findSpacer(t, p)
		}
	}
	fmt.Println("Now post-processing")
	s.postProcess(p) // remove non-unique entries from chromosome map
	x.postProcess(p) // remove non-unique entries from target map
	//if there is a target sequence, print out the target sequence
	if len(t) > 0 {
		s.findMatch(x, p)
		// no target sequence -> print out all unique target sites
	} else {
		s.printAll(p)
	}
	fmt.Println("Finished!!!")
}
