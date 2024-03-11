package GoDNA

import (
	"testing"
)

func TestContinuity(t *testing.T) {
	str := "GATCTATGTAAGGCCGGTTG"
	seq, _ := ToSeq(str)
	t.Log(Continuity(seq, 3))
}

func TestHairpin(t *testing.T) {
	str := "GTTCTATGTCCCACAGATCC"
	seq, _ := ToSeq(str)
	t.Log(Hairpin(seq, 6, 6, 3))
}

func TestHMeasure(t *testing.T) {
	str1 := "GATCTATGTAAGGCCGGTTG"
	str2 := "GTTCTATGTCCCACAGATCC"
	seq1, _ := ToSeq(str1)
	seq2, _ := ToSeq(str2)

	t.Log(HMeasure(seq2, seq2))
	t.Log(HMeasure(seq2, seq1))
}

func TestSimilarity(t *testing.T) {
	str1 := "GATCTATGTAAGGCCGGTTG"
	str2 := "GTTCTATGTCCCACAGATCC"
	seq1, _ := ToSeq(str1)
	seq2, _ := ToSeq(str2)

	t.Log(Similarity(seq1, seq2))
	t.Log(Similarity(seq2, seq1))
}

func TestMeltingTemperature(t *testing.T) {
	str1 := "CCTATCTCCTCATCCTCTAC"
	seq1, _ := ToSeq(str1)
	t.Log(MeltingTemperature(seq1))
}
