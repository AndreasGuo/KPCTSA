package DNAAlgorithm

import (
	"GoDNA/algorithm"
	"math"
)

type DNAPopulation struct {
	size         int
	varianceDim  int
	objectiveDim int
	lb, ub       float64
	individuals  []algorithm.Individual
	fit          FitFuncType
	zmin         []float64
}

func (pop *DNAPopulation) SetFitFunc(fit FitFuncType) {
	pop.fit = fit
}

func (pop *DNAPopulation) SetConfig(config *algorithm.Config) {
	pop.size = config.POPSIZE
	pop.varianceDim = config.DIM
	pop.objectiveDim = 5
	pop.lb = 0
	pop.ub = 4
}

func (pop *DNAPopulation) Size() int {
	return pop.size
}

func (pop *DNAPopulation) Init() {
	agents := make([]algorithm.Individual, pop.size)
	for i := range agents {
		agents[i] = CreateDNAAgent[float64](pop.varianceDim, pop.lb, pop.ub)
	}
}

func (pop *DNAPopulation) Append(invs []algorithm.Individual) {
	if invs != nil && len(invs) > 0 {
		for i := range invs {
			pop.individuals = append(pop.individuals, invs[i])
		}
	}
}

func (pop *DNAPopulation) Fit() [][]float64 {
	fits, ZMin := pop.fit(pop.individuals)
	pop.zmin = ZMin
	return fits
}

func (pop *DNAPopulation) ZMin() []float64 {
	return pop.zmin
}

func (pop *DNAPopulation) At(i int) algorithm.Individual {
	return pop.individuals[i]
}

func (pop *DNAPopulation) VarianceDim() int {
	return pop.varianceDim
}

func (pop *DNAPopulation) ObjectiveDim() int {
	return pop.objectiveDim
}

func (pop *DNAPopulation) UpdatePosition(i int, position []float64) {
	for i, v := range position {
		position[i] = math.Round(v)
	}
	pop.individuals[i].UpdatePosition(position)
}

func (pop *DNAPopulation) LB() float64 {
	return pop.lb
}

func (pop *DNAPopulation) UB() float64 {
	return pop.ub
}

func (pop *DNAPopulation) PostWork() {
	for i := range pop.Size() {
		pop.individuals[i].PostWork()
	}
}
