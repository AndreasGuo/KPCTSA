package algorithm

import (
	DNAType "GoDNA/algorithm/dnatype"
	"math"
	"math/rand"
)

// Grey Wolf Optimizer
type GWO struct {
	Pop          *DNAType.DNAPopulation
	MaxIteration int
}

func NewGWO(maxIteration int) Algorithm {
	return &GWO{Pop: nil, MaxIteration: maxIteration}
}

func (g *GWO) GetName() string {
	return "GWOX"
}

func (g *GWO) Initialize(pop *DNAType.DNAPopulation, inds ...*DNAType.DNAAgent) {
	pop.Init()
	if len(inds) > 0 {
		pop.Append(inds)
	}
	g.Pop = pop
}

func (g *GWO) Iteration(hyperPlaneNorm bool, origin bool, cd bool) *DNAType.DNAAgent {
	// logger := log.Default()
	// islog := false
	fits := g.Pop.Fit()
	ZMin := g.Pop.ZMin()
	bestIndex, argSorted := NDKPSort(fits, ZMin, g.Pop.Size(), hyperPlaneNorm, cd)

	bestIndividual := g.Pop.At(bestIndex)
	alpha := bestIndividual
	beta := g.Pop.At(argSorted[1])
	delta := g.Pop.At(argSorted[2])

	var l float64 = 0
	a := 2 - l*((2)/float64(g.MaxIteration))
	for it := 0; it < int(g.MaxIteration); it++ {
		oldPop := g.Pop.Clone()
		gamma := 0.5 - 0.4*(1-(float64(it)/float64(g.MaxIteration)))
		//extraIndex := oldPop.Size()
		for i := range oldPop.Size() {
			// if i == bestIndex {
			// 	continue
			// }
			x := oldPop.At(i).Variance()
			for j := range g.Pop.VarianceDim() {
				r1 := rand.Float64()
				r2 := rand.Float64()
				A1 := 2*a*r1 - a
				C1 := 2 * r2
				DAlpha := math.Abs(C1*float64(alpha.Seq[j]) - x[j])
				X1 := float64(alpha.Seq[j]) - A1*DAlpha

				r1 = rand.Float64()
				r2 = rand.Float64()
				A2 := 2*a*r1 - a
				C2 := 2 * r2
				DBeta := math.Abs(C2*float64(beta.Seq[j]) - x[j])
				X2 := float64(beta.Seq[j]) - A2*DBeta

				r1 = rand.Float64()
				r2 = rand.Float64()
				A3 := 2*a*r1 - a
				C3 := 2 * r2
				DDelta := math.Abs(C3*float64(delta.Seq[j]) - x[j])
				X3 := float64(beta.Seq[j]) - A3*DDelta

				if rand.Float64() < 0.1 {
					x[j] = (X1 + X2 + X3 + x[j] + C0(gamma)) / 4
				} else {
					x[j] = (X1 + X2 + X3) / 3
				}

				// boundary control
				x[j] = max(x[j], g.Pop.LB())
				x[j] = min(x[j], g.Pop.UB())
			}
			g.Pop.UpdatePosition(i, x)
		}
		l = l + 1

		g.Pop.PostWork()
		g.Pop.Join(oldPop)

		fits = g.Pop.Fit()
		ZMin = g.Pop.ZMin()
		bestIndex, selectedIndex := NDKPSort(fits, ZMin, oldPop.Size(), hyperPlaneNorm, cd)
		bestIndividual = g.Pop.At(bestIndex)
		g.Pop.Select(selectedIndex)
	}
	return bestIndividual
}
