package test

import (
	DNAType "GoDNA/algorithm/dnatype"
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
)

type PO struct {
	Pop          *DNAType.DNAPopulation
	MaxIteration int
}

func (po *PO) Initialize(pop *DNAType.DNAPopulation, inds ...*DNAType.DNAAgent) {
	pop.Init()
	if len(inds) > 0 {
		pop.Append(inds)
	}
	po.Pop = pop
}

func Levy(lh int) []float64 {
	beta := 1.5
	sigma := math.Gamma(1+beta) * math.Sin(math.Pi*beta/2)
	sigma /= math.Gamma(1+beta/2) * beta * math.Pow(1, (beta-1)/2)
	sigma = math.Pow(sigma, 1/beta)
	u := make([]float64, lh)
	for i := range u {
		u[i] = sigma * rand.NormFloat64() / math.Pow(math.Abs(rand.NormFloat64()), 1/beta)
	}
	return u
}

func Mean(pop *DNAType.DNAPopulation) []float64 {
	size := pop.Size()
	dim := pop.VarianceDim()
	out := make([]float64, dim)
	for i := range size {
		v := pop.At(i).Variance()
		for j := range dim {
			out[j] += v[j]
		}
	}

	for j := range dim {
		out[j] /= float64(size)
	}

	return out
}

func PrintDNASet(dnaSet []*DNAType.DNAAgent, fitChan *DNAType.FitChan) string {
	var result = ""
	for ind, inv := range dnaSet {
		DNAString, err := inv.String()
		if err != nil {
			panic("error while decoding")
		}
		result += strconv.Itoa(ind) + " " + DNAString

		// continuity
		fitChan.CtIn <- DNAType.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		ct := (<-fitChan.CtRe).Value
		result += " " + strconv.FormatFloat(ct, 'f', 4, 64)

		// hairpin
		fitChan.HpIn <- DNAType.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		hp := (<-fitChan.HpRe).Value
		result += " " + strconv.FormatFloat(hp, 'f', 4, 64)

		// hm
		hmList := make([]float64, len(dnaSet))
		for j, o := range dnaSet {
			fitChan.HmIn <- DNAType.SeqMapPair{Index1: ind, Index2: j, Seq1: inv.Represent(), Seq2: o.Represent()}
			hmList[j] = (<-fitChan.HmRe).Value
		}
		hm := sum(hmList)
		result += " " + strconv.FormatFloat(hm, 'f', 4, 64)

		// sm
		smList := make([]float64, len(dnaSet))
		for j, o := range dnaSet {
			if j != ind {
				fitChan.SmIn <- DNAType.SeqMapPair{Index1: ind, Index2: j, Seq1: inv.Represent(), Seq2: o.Represent()}
				smList[j] = (<-fitChan.SmRe).Value
			}
		}
		sm := sum(smList)
		result += " " + strconv.FormatFloat(sm, 'f', 4, 64)

		//mt
		fitChan.MtIn <- DNAType.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		mt := (<-fitChan.MtRe).Value
		result += " " + strconv.FormatFloat(mt, 'f', 4, 64) + "\n"
	}
	fmt.Println(result)
	return result
}

func sum[T int | float64](lt []T) T {
	var s T = 0
	for i := range lt {
		s += lt[i]
	}
	return s
}

func writeCSV(data [][]float64, behavior string, it int, origin bool, bestIdx int) {
	isOrigin := "origin"
	if !origin {
		isOrigin = "adjust"
	}
	filename := fmt.Sprintf("obj-%s-%d-%s.csv", behavior, it+1, isOrigin)
	file, err := os.OpenFile(filename, os.O_WRONLY|os.O_CREATE, 0666)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	defer writer.Flush()
	writer.Write([]string{"Continuity", "Hairpin", "HMeasure", "Similarity", "VMT", "Best"})
	for idx, value := range data {
		v := []string{}
		for _, iv := range value {
			v = append(v, strconv.FormatFloat(float64(iv), 'f', 2, 64))
		}
		if idx == bestIdx {
			v = append(v, "true")
		} else {
			v = append(v, "false")
		}
		writer.Write(v)
	}
}

func randomDNASet(dim, size int) []*DNAType.DNAAgent {
	dnaSet := []*DNAType.DNAAgent{}
	for i := 0; i < size; i++ {
		dnaSet = append(dnaSet, DNAType.CreateDNAAgent(dim, 0, 3))
	}
	return dnaSet
}

func agentFromStr(str string) *DNAType.DNAAgent {
	seq, _ := DNAType.ToSeq(str)
	agent := DNAType.DNAAgent{}
	// init raw variance
	agent.DNARowVariance = make([]float64, len(str))
	for i := 0; i < len(str); i++ {
		agent.DNARowVariance[i] = float64(seq[i])
	}
	agent.RepairAndToSeq()
	return &agent
}

func baseDNAs() []*DNAType.DNAAgent {
	dnas := []*DNAType.DNAAgent{}
	// CCTCCACTATCAGAACTACC
	// AACCTAACCTCCTCTCCATC
	agent1 := agentFromStr("CCTCCACTATCAGAACTACC")
	dnas = append(dnas, agent1)

	agent2 := agentFromStr("AACCTAACCTCCTCTCCATC")
	dnas = append(dnas, agent2)

	agent3 := DNAType.CreateDNAAgent(20, 0, 3)
	dnas = append(dnas, agent3)

	return dnas
}
