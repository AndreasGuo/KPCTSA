package algorithm

import (
	"GoDNA/DNAAnalysis"
	"fmt"
	"math/rand"
	"time"
)

type DNAAgent[T int | float64] struct {
	Variance []float64
	Seq      DNAAnalysis.Seq
	Ct       T
	Hp       T
	Hm       T
	Sm       T
	Mt       float64
}

func (dnaAgent *DNAAgent[T]) Objs() []float64 {
	objs := make([]float64, 5)
	objs[0] = float64(dnaAgent.Ct)
	objs[1] = float64(dnaAgent.Hp)
	objs[2] = float64(dnaAgent.Hm)
	objs[3] = float64(dnaAgent.Sm)
	objs[4] = float64(dnaAgent.Mt)
	return objs
}

func CreateDNAAgent[T int | float64](dim int, lb, ub float64) DNAAgent[T] {
	agent := DNAAgent[T]{}
	// init raw variance
	agent.Variance = make([]float64, dim)
	for i := 0; i < dim; i++ {
		agent.Variance[i] = rand.Float64()*(ub-lb) + lb
	}
	agent.RepairAndToSeq()
	return agent
}

func (dnaAgent *DNAAgent[T]) RepairAndToSeq() {
	// First control GC at 50%
	// and generate DNA seq.
	GCPosition := []int{} //make([]int, len(agent.variance))
	ATPosition := []int{} //make([]int, len(agent.variance))
	for i, value := range dnaAgent.Variance {
		if value == DNAAnalysis.C || value == DNAAnalysis.G {
			GCPosition = append(GCPosition, i)
		} else {
			ATPosition = append(ATPosition, i)
		}
	}
	if len(GCPosition) > len(dnaAgent.Variance)/2 {
		toRepairLen := len(GCPosition) - len(dnaAgent.Variance)/2
		shuffle(GCPosition)
		toRepair := GCPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				dnaAgent.Variance[toRepair[i]] = DNAAnalysis.T
			} else {
				dnaAgent.Variance[toRepair[i]] = DNAAnalysis.A
			}
		}
	}

	if len(ATPosition) > len(dnaAgent.Variance)/2 {
		toRepairLen := len(ATPosition) - len(dnaAgent.Variance)/2
		shuffle(ATPosition)
		toRepair := ATPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				dnaAgent.Variance[toRepair[i]] = DNAAnalysis.C
			} else {
				dnaAgent.Variance[toRepair[i]] = DNAAnalysis.G
			}
		}
	}

	for i, value := range dnaAgent.Variance {
		switch int(value) {
		case 0:
			dnaAgent.Seq[i] = DNAAnalysis.C
		case 1:
			dnaAgent.Seq[i] = DNAAnalysis.T
		case 2:
			dnaAgent.Seq[i] = DNAAnalysis.A
		case 3:
			dnaAgent.Seq[i] = DNAAnalysis.G
		default:
			panic(fmt.Sprintf("error while repair DNA seq, value %f do not defined!", value))
		}
	}
}

func shuffle(slice []int) {
	r := rand.New(rand.NewSource(time.Now().Unix()))
	for len(slice) > 0 {
		n := len(slice)
		randIndex := r.Intn(n)
		slice[n-1], slice[randIndex] = slice[randIndex], slice[n-1]
		slice = slice[:n-1]
	}
}

type DNAPopulation struct {
	individuals  []DNAAgent[float64]
	dim          int
	size         int
	maxIteration int
	lb, ub       float64
	//fitnessFunc  func([]Individual[T]) [][]float64
}

func (pop *DNAPopulation) Objs() [][]float64 {
	objs := make([][]float64, pop.size)
	for i := range objs {
		objs[i] = pop.individuals[i].Objs()
	}
	return objs
}

func InitDNAPopulation(dim, size, maxIteration int, lb, ub float64, agent DNAAgent[float64]) DNAPopulation {
	agents := make([]DNAAgent[float64], size)
	for i := range agents {
		agents[i] = CreateDNAAgent[float64](dim, lb, ub)
	}
	agents = append(agents, agent)
	pop := DNAPopulation{agents, dim, size + 1, maxIteration, lb, ub}
	return pop
}
