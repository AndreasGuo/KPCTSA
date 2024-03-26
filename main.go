package main

import (
	"GoDNA/DNAAnalysis"
	"GoDNA/algorithm"
	"fmt"
)

func main() {
	fitChan := DNAAnalysis.CreateWorker(100, 100, 10)
	defer fitChan.Close()

	var config = algorithm.DefaultConfig()
	var dnaSet = RandomDNASet(config.DIM, config.DNASIZE)

	for it := 0; it < config.DNASETITERATION; it++ {
		fmt.Print("DNASet iteration ", it+1, " ...")
		for index := range dnaSet {
			fitFunc := algorithm.FitnessCall(dnaSet, index, fitChan, config)
			alg := algorithm.PO{}
			alg.Initialize(config.DIM, config.POPSIZE, config.MAXIT,
				float64(config.LB), float64(config.UB), fitFunc, dnaSet[index])
			inv := alg.Iteration()

			dnaSet[index] = *inv
		}
		fmt.Println("Done")
		for ind, inv := range dnaSet {
			DNAString, err := inv.Seq.ToStr()
			if err != nil {
				panic("error while decoding")
			}
			fmt.Print(DNAString, " ")
			// continuity
			fitChan.CtIn <- DNAAnalysis.SeqMapSingle{ind, inv.Seq}
			fmt.Print((<-fitChan.CtRe).Value, " ")
			// hairpin
			fitChan.HpIn <- DNAAnalysis.SeqMapSingle{ind, inv.Seq}
			fmt.Print((<-fitChan.HpRe).Value, " ")
			// hm
			hmList := make([]float64, len(dnaSet))
			for j, o := range dnaSet {
				fitChan.HmIn <- DNAAnalysis.SeqMapPair{ind, j, inv.Seq, o.Seq}
				hmList[j] = (<-fitChan.HmRe).Value
			}
			fmt.Print(sum(hmList), " ")
			// sm
			smList := make([]float64, len(dnaSet))
			for j, o := range dnaSet {
				if j != ind {
					fitChan.SmIn <- DNAAnalysis.SeqMapPair{ind, j, inv.Seq, o.Seq}
					smList[j] = (<-fitChan.SmRe).Value
				}
			}
			fmt.Print(sum(smList), " ")
			//mt
			fitChan.MtIn <- DNAAnalysis.SeqMapSingle{ind, inv.Seq}
			fmt.Printf("%.4f\n", (<-fitChan.MtRe).Value)
		}
	}
}

func RandomDNASet(dim, size int) []algorithm.DNAAgent[float64] {
	dnaSet := []algorithm.DNAAgent[float64]{}
	for i := 0; i < size; i++ {
		dnaSet = append(dnaSet, algorithm.CreateDNAAgent[float64](dim, 0, 3))
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
