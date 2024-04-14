package DNAAlgorithm

import (
	"GoDNA/DNAAnalysis"
	"testing"
)

func TestDNAAgent_NoRunLength(t *testing.T) {
	str := "TCTCTCACACCTCTCTCTCT"
	agent := seqToAgent(str)
	t.Log(agent.Variance())
	agent.NoRunLength()
	t.Log(agent.Variance())
}

func seqToAgent(strSeq string) *DNAAgent {
	variance := make([]float64, len(strSeq))
	for i := range variance {
		switch strSeq[i] {
		case 'C':
			variance[i] = DNAAnalysis.C
		case 'T':
			variance[i] = DNAAnalysis.T
		case 'A':
			variance[i] = DNAAnalysis.A
		case 'G':
			variance[i] = DNAAnalysis.G
		default:
			panic("啥玩意啊" + string(strSeq[i]))
		}
	}
	return &DNAAgent{DNARowVariance: variance}
}
