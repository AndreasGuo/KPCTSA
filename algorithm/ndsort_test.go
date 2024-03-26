package algorithm

import "testing"

func TestNDSort(t *testing.T) {
	data := [][]int{
		{39, 27, 11, 42},
		{10, 93, 91, 90},
		{54, 78, 56, 89},
		{24, 64, 20, 65},
	}
	t.Log(NDSort(data))
}
