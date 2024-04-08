package DNAAlgorithm

import (
	"GoDNA/DNAAnalysis"
	"fmt"
	"math"
	"math/rand"
	"time"
)

type DNAAgent struct {
	DNARowVariance []float64
	Seq            DNAAnalysis.Seq
	Ct             float64
	Hp             float64
	Hm             float64
	Sm             float64
	Mt             float64
}

func (dnaAgent *DNAAgent) String() (string, error) {
	return dnaAgent.Seq.ToStr()
}

func (dnaAgent *DNAAgent) Represent() []int {
	return dnaAgent.Seq
}

func (dnaAgent *DNAAgent) UpdatePosition(position []float64) {
	dnaAgent.DNARowVariance = position
}

func (dnaAgent *DNAAgent) Objs() []float64 {
	objs := make([]float64, 5)
	objs[0] = float64(dnaAgent.Ct)
	objs[1] = float64(dnaAgent.Hp)
	objs[2] = float64(dnaAgent.Hm)
	objs[3] = float64(dnaAgent.Sm)
	objs[4] = float64(dnaAgent.Mt)
	return objs
}

func (dnaAgent *DNAAgent) Variance() []float64 {
	return dnaAgent.DNARowVariance
}

func CreateDNAAgent(dim int, lb, ub float64) *DNAAgent {
	agent := DNAAgent{}
	// init raw variance
	agent.DNARowVariance = make([]float64, dim)
	for i := 0; i < dim; i++ {
		agent.DNARowVariance[i] = math.Round(rand.Float64()*(ub-lb) + lb)
	}
	agent.RepairAndToSeq()
	return &agent
}

func (dnaAgent *DNAAgent) PostWork() {
	dnaAgent.RepairAndToSeq()
}

func (dnaAgent *DNAAgent) RepairAndToSeq() {
	dnaAgent.fixGCContent()
	dnaAgent.runlongConstraint()
	// generate DNA seq.
	dnaAgent.Seq = convSeq(dnaAgent.DNARowVariance)
}

func (dnaAgent *DNAAgent) fixGCContent() {
	// First control GC at 50%
	GCPosition := []int{} //make([]int, len(agent.variance))
	ATPosition := []int{} //make([]int, len(agent.variance))
	for i, value := range dnaAgent.DNARowVariance {
		if value == DNAAnalysis.C || value == DNAAnalysis.G {
			GCPosition = append(GCPosition, i)
		} else {
			ATPosition = append(ATPosition, i)
		}
	}
	if len(GCPosition) > len(dnaAgent.DNARowVariance)/2 {
		toRepairLen := len(GCPosition) - len(dnaAgent.DNARowVariance)/2
		shuffle(GCPosition)
		toRepair := GCPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				dnaAgent.DNARowVariance[toRepair[i]] = DNAAnalysis.T
			} else {
				dnaAgent.DNARowVariance[toRepair[i]] = DNAAnalysis.A
			}
		}
	}

	if len(ATPosition) > len(dnaAgent.DNARowVariance)/2 {
		toRepairLen := len(ATPosition) - len(dnaAgent.DNARowVariance)/2
		shuffle(ATPosition)
		toRepair := ATPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				dnaAgent.DNARowVariance[toRepair[i]] = DNAAnalysis.C
			} else {
				dnaAgent.DNARowVariance[toRepair[i]] = DNAAnalysis.G
			}
		}
	}
}

func (DNAAgent *DNAAgent) runlongConstraint() {
	// Assure no continuous identical bases
	// TODO
}

func convSeq(variance []float64) DNAAnalysis.Seq {
	// convert float variance array to DNA Seq
	seq := make(DNAAnalysis.Seq, len(variance))
	for i, value := range variance {
		switch int(value) {
		case 0:
			seq[i] = DNAAnalysis.C
		case 1:
			seq[i] = DNAAnalysis.T
		case 2:
			seq[i] = DNAAnalysis.A
		case 3:
			seq[i] = DNAAnalysis.G
		default:
			panic(fmt.Sprintf("error while repair DNA seq, value %f do not defined!", value))
		}
	}
	return seq
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
