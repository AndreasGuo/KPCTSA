package main

import (
	"GoDNA/DNAAlgorithm"
	"GoDNA/DNAAnalysis"
	"GoDNA/algorithm"
	"fmt"
	"os"
	"strconv"
)

func main() {
	fitChan := DNAAnalysis.CreateWorker(100, 100, 10)
	defer fitChan.Close()

	var config = DNAAlgorithm.DefaultConfig()
	var dnaSet = RandomDNASet(config.DIM, config.DNASIZE)

	result := ""
	printDNASet(dnaSet, fitChan, result)

	for it := 0; it < config.DNASETITERATION; it++ {
		fmt.Println("DNASet iteration ", it+1, "/", config.DNASETITERATION)

		// init objs
		for i := range dnaSet {
			fitFunc := DNAAlgorithm.FitnessCall(dnaSet, i, fitChan, config)
			singleInvSlice := []algorithm.Individual{dnaSet[i]}
			fits, _ := fitFunc(singleInvSlice)
			dnaSet[i].SetObjs(fits[0])
		}

		index := chooseInvToOpt(dnaSet)
		fmt.Println("To Optimize: ", index)

		//for index := range dnaSet {
		fitFunc := DNAAlgorithm.FitnessCall(dnaSet, index, fitChan, config)
		alg := algorithm.PO{Pop: nil, MaxIteration: config.MAXIT}
		// alg := algorithm.BKA{Pop: nil, MaxIteration: config.MAXIT}
		pop := new(DNAAlgorithm.DNAPopulation)
		pop.SetConfig(config)
		pop.SetFitFunc(fitFunc)
		alg.Initialize(pop, dnaSet[index])
		inv := alg.Iteration()
		dnaSet[index] = inv
		//}

		fmt.Println("\rDone")
		result = printDNASet(dnaSet, fitChan, result)

	}
	fileName := "result.txt"
	os.WriteFile(fileName, []byte(result), 0644)
}

func printDNASet(dnaSet []algorithm.Individual, fitChan *DNAAnalysis.FitChan, result string) string {
	result = ""
	for ind, inv := range dnaSet {
		DNAString, err := inv.String()
		if err != nil {
			panic("error while decoding")
		}
		result += strconv.Itoa(ind) + " " + DNAString
		//fmt.Print(DNAString, " ")
		// continuity
		fitChan.CtIn <- DNAAnalysis.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		ct := (<-fitChan.CtRe).Value
		result += " " + strconv.FormatFloat(ct, 'f', 4, 64)
		//fmt.Print((<-fitChan.CtRe).Value, " ")
		// hairpin
		fitChan.HpIn <- DNAAnalysis.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		hp := (<-fitChan.HpRe).Value
		result += " " + strconv.FormatFloat(hp, 'f', 4, 64)
		//fmt.Print((<-fitChan.HpRe).Value, " ")
		// hm
		hmList := make([]float64, len(dnaSet))
		for j, o := range dnaSet {
			fitChan.HmIn <- DNAAnalysis.SeqMapPair{Index1: ind, Index2: j, Seq1: inv.Represent(), Seq2: o.Represent()}
			hmList[j] = (<-fitChan.HmRe).Value
		}
		//fmt.Print(sum(hmList), " ")
		hm := sum(hmList)
		result += " " + strconv.FormatFloat(hm, 'f', 4, 64)
		// sm
		smList := make([]float64, len(dnaSet))
		for j, o := range dnaSet {
			if j != ind {
				fitChan.SmIn <- DNAAnalysis.SeqMapPair{Index1: ind, Index2: j, Seq1: inv.Represent(), Seq2: o.Represent()}
				smList[j] = (<-fitChan.SmRe).Value
			}
		}
		//fmt.Print(sum(smList), " ")
		sm := sum(smList)
		result += " " + strconv.FormatFloat(sm, 'f', 4, 64)
		//mt
		fitChan.MtIn <- DNAAnalysis.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		//fmt.Printf("%.4f\n", (<-fitChan.MtRe).Value)
		mt := (<-fitChan.MtRe).Value
		result += " " + strconv.FormatFloat(mt, 'f', 4, 64) + "\n"
		// dnaSet[ind].SetObjs([]float64{ct, hp, hm, sm, mt})
	}
	fmt.Println()
	fmt.Println(result)
	return result
}

func RandomDNASet(dim, size int) []algorithm.Individual {
	dnaSet := []algorithm.Individual{}
	for i := 0; i < size; i++ {
		dnaSet = append(dnaSet, DNAAlgorithm.CreateDNAAgent(dim, 0, 3))
	}
	return dnaSet
}

func sum[T int | float64](lt []T) T {
	var s T = 0
	for i := range lt {
		s += lt[i]
	}
	return s
}

func chooseInvToOpt(dnaSet []algorithm.Individual) (idx int) {
	if len(dnaSet) == 0 {
		panic("no inv")
	}
	ZMin := make([]float64, len(dnaSet[0].Objs()))
	for j := range ZMin {
		ZMin[j] = 400
	}
	for _, dna := range dnaSet {
		objs := dna.Objs()
		for i := range len(ZMin) {
			ZMin[i] = min(ZMin[i], objs[i])
		}
	}
	distance := make([]float64, len(dnaSet))
	for i, dna := range dnaSet {
		objs := dna.Objs()
		distance[i] = 1
		for j := range ZMin {
			//distance[i] += math.Pow(objs[j]-ZMin[j], 2)
			distance[i] *= max(1, objs[j]-ZMin[j])
		}
	}

	dix := 0.0
	for i := range distance {
		if distance[i] > dix {
			idx = i
			dix = distance[i]
		}
	}
	return idx
}
