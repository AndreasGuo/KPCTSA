package algorithm

// NDSort
// Efficient non-dominated sort with sequential search (ENS-SS)
// @Input objs matrix of object values, size*object
// @Output front no of each individual
func NDSort[T int | float64](objs [][]T) []int {
	var maxFrontNo = -1
	var frontNo = make([]int, len(objs))
	for i := range frontNo {
		frontNo[i] = -1
	}

	for sorted(frontNo) < len(frontNo) {
		maxFrontNo += 1
		marked := false
		for i := range objs {
			if frontNo[i] == -1 {
				dominated := false

				for j := 0; j < len(objs); j++ {
					if i == j {
						continue
					}
					if !marked && frontNo[j] == -1 && dominate(objs[j], objs[i]) {
						dominated = true
						break
					} else if (frontNo[j] == maxFrontNo) && dominate(objs[j], objs[i]) {
						dominated = true
						break
					}
				}

				if !dominated {
					frontNo[i] = maxFrontNo
					marked = true
				}
			}
		}
	}

	return frontNo
}

func sorted(frontNo []int) int {
	count := 0
	for _, v := range frontNo {
		if v != -1 {
			count += 1
		}
	}
	return count
}

func dominate[T int | float64](lhs, rhs []T) bool {
	count := 0
	atLeastLess := false
	for i := range lhs {
		if lhs[i] <= rhs[i] {
			count += 1
		}
		if lhs[i] < rhs[i] {
			atLeastLess = true
		}
	}
	return (count == len(lhs)) && atLeastLess
}
