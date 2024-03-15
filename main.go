package main

import "fmt"

const DIM = 20
const DNASIZE = 7
const POPSIZE = 20
const MAXIT = 100
const LB = 0
const UB = 3
const DNASETITERATION = 1000

func main() {
	fitChan := CreateWorker(100, 100, 10)
	defer fitChan.Close()

	var dnaSet = RandomDNASet(DIM, DNASIZE)

	for it := 0; it < DNASETITERATION; it++ {
		fmt.Print("DNASet iteration ", it+1, " ...")
		for index := range dnaSet {
			pop := CreatePopulation(DIM, POPSIZE, MAXIT, LB, UB)
			fitFunc := FitnessCall(dnaSet, index, fitChan)
			inv := pop.UpdatePopulation(fitFunc)
			if it == 0 || inv.fitness < dnaSet[index].fitness {
				dnaSet[index] = inv
			}
		}
		fmt.Println("Done")
		for _, inv := range dnaSet {
			DNAString, err := inv.seq.ToStr()
			if err != nil {
				panic("error while decoding")
			}
			fmt.Println(DNAString)
		}
	}

}

func RandomDNASet(dim, size uint) []individual {
	pop := CreatePopulation(dim, size, 0, 0, 3)
	//pop.Fitness(fitFunc)
	return pop.individuals
}
