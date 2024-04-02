package algorithm

import (
	"errors"

	"gonum.org/v1/gonum/mat"
)

const w float64 = 1e6

func bestIndex(points [][]float64, ZMin []float64) (int, error) {
	n := len(points)
	if n == 0 {
		return 0, errors.New("There is no points")
	}
	m := len(points[0])
	if m == 0 {
		return 0, errors.New("Dimension is zero")
	}

	extrems := []float64{}

	localPoints := make([][]float64, n)
	copy(localPoints, points)

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			localPoints[i][j] -= ZMin[j]
		}
	}

	// ASF
	for j := 0; j < m; j++ {
		minASFIndexDimM := 0
		minASF := -1.0
		for i := 0; i < n; i++ {
			maxASF := float64(0)
			for jj := 0; jj < m; jj++ {
				if j == jj {
					maxASF = max(maxASF, localPoints[i][jj])
				} else {
					maxASF = max(maxASF, localPoints[i][jj]*w)
				}

			}
			if minASFIndexDimM == 0 || maxASF < minASF {
				minASF = maxASF
				minASFIndexDimM = i
			}
		}
		extrems = append(extrems, localPoints[minASFIndexDimM]...)

	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			localPoints[i][j] -= ZMin[j]
		}
	}

	one := make([]float64, m)
	for i := range one {
		one[i] = 1
	}
	ones := mat.NewDense(m, 1, one)

	extMtrix := mat.NewDense(m, m, extrems)
	var aInv mat.Dense
	aInv.Solve(extMtrix, ones)

	intercept := make([]float64, m)
	for j := 0; j < m; j++ {
		intercept[j] = aInv.At(0, j)
	}
	for i := range intercept {
		if intercept[i] == 0 {
			for j := range localPoints {
				if intercept[i] < localPoints[j][i] {
					intercept[i] = localPoints[j][i]
				}
			}
		} else {
			intercept[i] = 1.0 / intercept[i]
		}
	}

	// 点到平面距离，并不是真的距离，只是不影响比较大小
	minDistance := -1.1
	index := 0
	for i := range localPoints {
		distance := float64(0)
		for j := 0; j < m; j++ {
			distance += intercept[j] * localPoints[i][j]
		}
		if index == 0 || distance < minDistance {
			index = i
			minDistance = distance
		}
	}
	return index, nil
}
