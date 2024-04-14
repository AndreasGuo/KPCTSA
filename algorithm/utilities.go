package algorithm

func sum[T int | float64](in []T) T {
	var s T = 0
	for i := range in {
		s += in[i]
	}
	return s
}
