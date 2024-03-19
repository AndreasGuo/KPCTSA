package main

import "fmt"

const DIM = 20
const DNASIZE = 7
const POPSIZE = 30
const MAXIT = 500
const LB = 0
const UB = 3
const DNASETITERATION = 200
const MINVALUE = 2e-2

func main() {
	fitChan := CreateWorker(100, 100, 10)
	defer fitChan.Close()

	var dnaSet = RandomDNASet(DIM, DNASIZE)

	for it := 0; it < DNASETITERATION; it++ {
		fmt.Print("DNASet iteration ", it+1, " ...")
		for index := range dnaSet {
			pop := CreatePopulation(DIM, POPSIZE, MAXIT, LB, UB, &(dnaSet[index]))
			fitFunc := FitnessCall(dnaSet, index, fitChan)
			inv := pop.UpdatePopulation(fitFunc)

			dnaSet[index] = inv
		}
		fmt.Println("Done")
		for ind, inv := range dnaSet {
			DNAString, err := inv.seq.ToStr()
			if err != nil {
				panic("error while decoding")
			}
			fmt.Print(DNAString, " ")
			// continuity
			fitChan.ctIn <- seqMapSingle{ind, inv.seq}
			fmt.Print((<-fitChan.ctRe).value, " ")
			// hairpin
			fitChan.hpIn <- seqMapSingle{ind, inv.seq}
			fmt.Print((<-fitChan.hpRe).value, " ")
			// hm
			hmList := make([]float64, len(dnaSet))
			for j, o := range dnaSet {
				fitChan.hmIn <- seqMapPair{ind, j, inv.seq, o.seq}
				hmList[j] = (<-fitChan.hmRe).value
			}
			fmt.Print(sum(hmList), " ")
			// sm
			smList := make([]float64, len(dnaSet))
			for j, o := range dnaSet {
				if j != ind {
					fitChan.smIn <- seqMapPair{ind, j, inv.seq, o.seq}
					smList[j] = (<-fitChan.smRe).value
				}
			}
			fmt.Print(sum(smList), " ")
			//mt
			fitChan.mtIn <- seqMapSingle{ind, inv.seq}
			fmt.Printf("%.4f\n", (<-fitChan.mtRe).value)
		}
	}

}

func RandomDNASet(dim, size uint) []individual {
	pop := CreatePopulation(dim, size, 0, 0, 3, nil)
	//pop.Fitness(fitFunc)
	return pop.individuals
}
