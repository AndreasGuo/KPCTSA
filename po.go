package main

import (
	"GoDNA/DNAAnalysis"
	"math"
	"math/rand"
	"time"
)

type individual struct {
	// Since variance is store in type float64, it is necessary to define another var to store DNA base type.
	// Type cannot be converted inherent.
	variance []float64
	seq      DNAAnalysis.Seq
	fitness  float64
}

type population struct {
	individuals  []individual
	dim          uint
	size         uint
	maxIteration uint
	lb, ub       float64
}

func CreatePopulation(dim, size uint, maxIteration uint, lb, ub float64) *population {
	pop := new(population)
	pop.dim = dim
	pop.size = size
	pop.lb, pop.ub = lb, ub
	pop.maxIteration = maxIteration

	pop.individuals = make([]individual, size)

	return pop
}

func levy[T int | float64](lh T) float64 {
	w := 1 + (float64(lh)-1)/4
	l := math.Pow(math.Sin(math.Pi*w), 2)
	// A step is ignored here which PO only use 1 dimension version of levy
	l += math.Pow(float64(lh-1), 2) * (1 + math.Pow(math.Sin(2*math.Pi*float64(lh)), 2))

	return l
}

func (pop *population) Fitness(fitnessFunc func([]float64) float64) int {
	bestIndex := -1
	for ind, inv := range pop.individuals {
		inv.fitness = fitnessFunc(inv.variance)
		if bestIndex == -1 || inv.fitness < pop.individuals[bestIndex].fitness {
			bestIndex = ind
		}
	}
	return bestIndex
}

func mean[T int | float64](x []T) float64 {
	var sum T = 0
	for i := 0; i < len(x); i++ {
		sum += x[i]
	}
	return float64(sum) / float64(len(x))
}

func st0(x, gbest []float64, dim, it, maxIt int) {
	ldim := levy(dim)
	meanx := mean(x)
	r := rand.Float64()
	for i := 0; i < dim; i++ {
		// X_new(j, :) = #inmatlab //(X(j, :) - GBestX) .* Levy(dim) + rand(1) * mean(X(j, :)) * (1 - i / Max_iter) ^ (2 * i / Max_iter);
		x[i] = (x[i]-gbest[i])*ldim + r*meanx*math.Pow(1-float64(it)/float64(maxIt), 2*float64(it)/float64(maxIt))
	}
}

func st1(x, gbest []float64, dim, it, maxIt int) {
	// X_new(j, :) = X(j, :) + GBestX .* Levy(dim) + randn() * (1 - i / Max_iter) * ones(1, dim);
	r := rand.NormFloat64()
	ldim := levy(dim)
	for i := 0; i < dim; i++ {
		x[i] = x[i] + gbest[i]*ldim + r*(1-float64(it)/float64(maxIt))
	}
}

func st2(x, gbest []float64, dim, it, maxIt int, alpha float64) {
	/*H = rand(1);
	if H < 0.5
		X_new(j, :) = X(j, :) + alpha * (1 - i / Max_iter) * (X(j, :) - mean(X(j, :)));
	else
	X_new(j, :) = X(j, :) + alpha * (1 - i / Max_iter) * exp(-j / (rand(1) * Max_iter));
	end*/
	h := rand.Float64()
	meanx := mean(x)
	if h < 0.5 {
		for i := 0; i < dim; i++ {
			x[i] += alpha * (1 - float64(it)/float64(maxIt)) * (x[i] - meanx)
		}
	} else {
		for i := 0; i < dim; i++ {
			x[i] += alpha * (1 - float64(it)/float64(maxIt)) * math.Exp(float64(-i)/(rand.Float64()*float64(maxIt)))
		}
	}
}

// X_new(j, :) = X(j, :) + rand() * cos((pi *i )/ (2 * Max_iter)) * (GBestX - X(j, :)) - cos(sita) * (i / Max_iter) ^ (2 / Max_iter) * (X(j, :) - GBestX);
func st3(x, gbest []float64, dim, it, maxIt int, sita float64) {
	r := rand.Float64()
	for i := 0; i < dim; i++ {
		x[i] += r*math.Cos(math.Pi*float64(it))/(2*float64(maxIt))*(gbest[i]-x[i]) -
			math.Cos(sita)*math.Pow((float64(it)/float64(maxIt)), 2/float64(maxIt))*
				x[i] - gbest[i]
	}
}

func repair(inv individual) {
	// First control GC at 50%
	// and generate DNA seq.
	GCPosition := make([]int, len(inv.variance))
	ATPosition := make([]int, len(inv.variance))
	for i, value := range inv.variance {
		if value == DNAAnalysis.C || value == DNAAnalysis.G {
			GCPosition = append(GCPosition, i)
		} else {
			ATPosition = append(ATPosition, i)
		}
	}
	if len(GCPosition) > len(inv.variance)/2 {
		toRepairLen := len(GCPosition) - len(inv.variance)/2
		shuffle(GCPosition)
		toRepair := GCPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				inv.variance[toRepair[i]] = DNAAnalysis.T
			} else {
				inv.variance[toRepair[i]] = DNAAnalysis.A
			}
		}
	}

	if len(ATPosition) > len(inv.variance)/2 {
		toRepairLen := len(ATPosition) - len(inv.variance)/2
		shuffle(ATPosition)
		toRepair := ATPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				inv.variance[toRepair[i]] = DNAAnalysis.C
			} else {
				inv.variance[toRepair[i]] = DNAAnalysis.G
			}
		}
	}

	for i, value := range inv.variance {
		switch value {
		case 0:
			inv.seq[i] = DNAAnalysis.C
		case 1:
			inv.seq[i] = DNAAnalysis.T
		case 2:
			inv.seq[i] = DNAAnalysis.A
		case 3:
			inv.seq[i] = DNAAnalysis.G
		default:
			panic("error while repair DNA seq, value do not defined!")
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

func (pop *population) UpdatePopulation(fitnessFunc func([]float64) float64) individual {
	bestIndex := pop.Fitness(fitnessFunc)
	gbest := pop.individuals[bestIndex].variance
	alpha := rand.Float64() / 5
	sita := rand.Float64() * math.Pi

	for it := 0; it < int(pop.maxIteration); it++ {
		for i := 0; i < int(pop.size); i++ {
			st := rand.Intn(4)
			switch st {
			case 0:
				st0(pop.individuals[i].variance, gbest, int(pop.dim), it, int(pop.maxIteration))
			case 1:
				st1(pop.individuals[i].variance, gbest, int(pop.dim), it, int(pop.maxIteration))
			case 2:
				st2(pop.individuals[i].variance, gbest, int(pop.dim), it, int(pop.maxIteration), alpha)
			case 3:
				st3(pop.individuals[i].variance, gbest, int(pop.dim), it, int(pop.maxIteration), sita)
			}
			//round and boundary control
			for j := 0; j < int(pop.dim); j++ {
				pop.individuals[i].variance[j] = math.Round(pop.individuals[i].variance[j])
				pop.individuals[i].variance[j] = max(pop.individuals[i].variance[j], pop.lb)
				pop.individuals[i].variance[j] = min(pop.individuals[i].variance[j], pop.ub)
			}
			repair(pop.individuals[i])
		}
		bestIndex = pop.Fitness(fitnessFunc)
	}
	return pop.individuals[bestIndex]
}
