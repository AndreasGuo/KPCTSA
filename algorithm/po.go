package algorithm

import (
	"fmt"
	"math"
	"math/rand"
)

type PO struct {
	Pop          Population
	MaxIteration int
}

func (po *PO) Initialize(pop Population, inds ...Individual) {
	pop.Init()
	if inds != nil && len(inds) > 0 {
		pop.Append(inds)
	}
	po.Pop = pop
}

// PO + NDSort + Knee Point
func (po *PO) Iteration() Individual {
	fits := po.Pop.Fit()
	ZMin := po.Pop.ZMin()
	selectedIndex := NDKPSort(fits, ZMin)

	gbest := po.Pop.At(selectedIndex).Variance()
	alpha := rand.Float64() / 5
	sita := rand.Float64() * math.Pi

	for it := 0; it < int(po.MaxIteration); it++ {
		fmt.Print("\r", it)
		for i := 0; i < int(po.Pop.Size()); i++ {
			st := rand.Intn(4)
			x := po.Pop.At(i).Variance()
			switch st {
			case 0:
				st0(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration)
			case 1:
				st1(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration)
			case 2:
				st2(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration, alpha)
			case 3:
				st3(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration, sita)
			}

			//boundary control
			for j := 0; j < po.Pop.VarianceDim(); j++ {
				x[j] = max(x[j], po.Pop.LB())
				x[j] = min(x[j], po.Pop.UB())
			}
			po.Pop.UpdatePosition(i, x)
			//po.pop.individuals[i].RepairAndToSeq()
		}
		po.Pop.PostWork()
		fits = po.Pop.Fit()
		ZMin = po.Pop.ZMin()
		selectedIndex = NDKPSort(fits, ZMin)
	}

	return po.Pop.At(selectedIndex)
}

func levy(lh int) []float64 {
	beta := 1.5
	sigma := math.Gamma(1+beta) * math.Sin(math.Pi*beta/2)
	sigma /= math.Gamma(1+beta/2) * beta * math.Pow(1, (beta-1)/2)
	sigma = math.Pow(sigma, 1/beta)
	u := make([]float64, lh)
	for i := range u {
		u[i] = rand.NormFloat64() / math.Pow(math.Abs(rand.NormFloat64()), 1/beta)
	}
	return u
}
func st0(x, gbest []float64, dim, it, maxIt int) {
	levyDim := levy(dim)
	meanx := mean(x)
	r := rand.Float64()
	for i := 0; i < dim; i++ {
		// X_new(j, :) = #inmatlab //(X(j, :) - GBestX) .* Levy(dim) + rand(1) * mean(X(j, :)) * (1 - i / Max_iter) ^ (2 * i / Max_iter);
		x[i] = (x[i]-gbest[i])*levyDim[i] + r*meanx*math.Pow(1-float64(it)/float64(maxIt), 2*float64(it)/float64(maxIt))
	}
}

func st1(x, gbest []float64, dim, it, maxIt int) {
	// X_new(j, :) = X(j, :) + GBestX .* Levy(dim) + randn() * (1 - i / Max_iter) * ones(1, dim);
	r := rand.NormFloat64()
	levyDim := levy(dim)
	for i := 0; i < dim; i++ {
		x[i] = x[i] + gbest[i]*levyDim[i] + r*(1-float64(it)/float64(maxIt))
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

func NDKPSort(fits [][]float64, ZMin []float64) int {
	NDFront := NDSort(fits)
	zeroPosition := findPosition(NDFront, 0)
	zeroFront := extractMatrix(fits, zeroPosition)
	bestIndex, _ := bestIndex(zeroFront, ZMin)
	return zeroPosition[bestIndex]
}

func mean[T int | float64](x []T) float64 {
	var sum T = 0
	for i := 0; i < len(x); i++ {
		sum += x[i]
	}
	return float64(sum) / float64(len(x))
}

func findPosition(lhs []int, rhs int) []int {
	positions := []int{}
	for i := range lhs {
		if lhs[i] == rhs {
			positions = append(positions, i)
		}
	}
	return positions
}

func extractMatrix(fits [][]float64, positions []int) [][]float64 {
	out := [][]float64{}
	for i := range positions {
		out = append(out, fits[positions[i]])
	}
	return out
}
