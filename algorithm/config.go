package algorithm

type Config struct {
	DIM             int
	DNASIZE         int
	POPSIZE         int
	MAXIT           int
	LB              int
	UB              int
	DNASETITERATION int
	MINVALUE        float64
}

func DefaultConfig() *Config {
	config := Config{20, 7, 30, 500, 0, 3, 200, 2e-2}
	return &config
}
