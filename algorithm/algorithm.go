package algorithm

type Algorithm[T int | float64] interface {
	// Initialize a population, with fitness calculated
	// The fitness function, it takes all individual as input and calculate fitness parallel
	// and returns the size*objs matrix of fitness values.
	Initialize(dim, size, maxIteration int, lb, ub float64, fitnessFunc func([]Individual[float64]) [][]float64)
	Iteration() Individual[T]
}

type Individual[T int | float64] interface {
	// Objs Returns object list
	Objs() []float64
}

type Population[T int | float64] interface {
	// Objs Returns object matrix
	Objs() [][]float64
}
