package algorithm

type anyFuncType func(...interface{}) any

type Algorithm[T int | float64] interface {
	// Initialize a population, with fitness calculated
	// The fitness function, it takes all individual as input and calculate fitness parallel
	// and returns the size*objs matrix of fitness values.
	Initialize(pop Population, inds ...Individual)
	Iteration() Individual
}

//type Population[T int | float64] struct {
//	Size        int
//	Lb, Ub      T
//	VarianceDim int
//	ObjectDim   int
//	Pop         []Individual
//}

type Population interface {
	Size() int
	Init()
	Append([]Individual)
	Fit() [][]float64
	ZMin() []float64
	At(int) Individual
	VarianceDim() int
	ObjectiveDim() int
	UpdatePosition(int, []float64)
	LB() float64
	UB() float64
	PostWork()
}

type Individual interface {
	// Objs Returns object list
	Objs() []float64
	Variance() []float64
	UpdatePosition([]float64)
	PostWork()
	Represent() []int
	String() (string, error)
}
