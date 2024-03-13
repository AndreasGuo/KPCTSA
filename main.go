package main

import (
	"GoDNA/DNAAnalysis"
	"fmt"
)

func main() {
	str1 := "GATCTATGTAAGGCCGGTTG"
	str2 := "GTTCTATGTCCCACAGATCC"
	seq1, _ := DNAAnalysis.ToSeq(str1)
	seq2, _ := DNAAnalysis.ToSeq(str2)

	continuityChan, continuityResult, hairpinChan, hairpinResult, hmChan, hmResult, smChan, smResult, mtChan, mtResult := CreateWorker(100, 100, 10)
	defer close(continuityChan)
	defer close(continuityResult)
	defer close(hairpinChan)
	defer close(hairpinResult)
	defer close(hmChan)
	defer close(hmResult)
	defer close(smChan)
	defer close(smResult)
	defer close(mtChan)
	defer close(mtResult)

	var in1 = seqMapSingle{0, seq1}
	var in2 = seqMapSingle{1, seq2}
	hairpinChan <- in1
	fmt.Println(<-hairpinResult)
	hairpinChan <- in2
	fmt.Println(<-hairpinResult)

	var pair = seqMapPair{0, 1, seq1, seq2}
	hmChan <- pair
	fmt.Println(<-hmResult)
	smChan <- pair
	fmt.Println(<-smResult)

	mtChan <- in1
	fmt.Println(<-mtResult)
	mtChan <- in2
	fmt.Println(<-mtResult)
}
