package main

import (
	"GoDNA/DNAAnalysis"
	"math"
	"slices"
	"sync"
)

var minCt, minHp, minHm, minSm = 400.0, 400.0, 400.0, 400.0

func FitnessCall(dnaSet DNASet, index int, fitChan *FitChan) func(individuals []individual) int {
	seqSet := make([]DNAAnalysis.Seq, len(dnaSet))
	for i := range seqSet {
		seqSet[i] = make(DNAAnalysis.Seq, DIM)
	}
	for i := range dnaSet {
		copy(seqSet[i], dnaSet[i].seq)
	}

	var mtValues = make([]float64, DNASIZE)
	for i, seq := range seqSet {
		fitChan.mtIn <- seqMapSingle{i, seq}
		re := <-fitChan.mtRe
		mtValues[re.index] = re.value
	}

	return func(invs []individual) int {
		go func() {
			for i := range invs {
				fitChan.ctIn <- seqMapSingle{i, invs[i].seq}
			}
		}()
		go func() {
			for i := range invs {
				fitChan.hpIn <- seqMapSingle{i, invs[i].seq}
			}
		}()
		go func() {
			for i := range invs {
				for j := range seqSet {
					if j == index {
						fitChan.hmIn <- seqMapPair{i, j, invs[i].seq, invs[i].seq}
					} else {
						// 正常的算法
						//fitChan.hmIn <- seqMapPair{i, j, invs[i].seq, seqSet[j]}
						// 交换前后
						fitChan.hmIn <- seqMapPair{i, j, seqSet[j], invs[i].seq}
					}
				}
			}
		}()
		go func() {
			for i := range invs {
				for j := range seqSet {
					if j == index {
						continue
					} else {
						// 正常的算法
						//fitChan.smIn <- seqMapPair{i, j, invs[i].seq, seqSet[j]}
						// 交换前后
						fitChan.smIn <- seqMapPair{i, j, seqSet[j], invs[i].seq}
					}
				}
			}
		}()
		go func() {
			for i := range invs {
				fitChan.mtIn <- seqMapSingle{i, invs[i].seq}
			}
		}()

		continuityList := make([]float64, len(invs))
		hairpinList := make([]float64, len(invs))
		hmTable := make([][]float64, len(invs))
		for i := range hmTable {
			hmTable[i] = make([]float64, len(dnaSet))
		}
		smTable := make([][]float64, len(invs))
		for i := range smTable {
			smTable[i] = make([]float64, len(dnaSet))
		}
		mtList := make([]float64, len(invs))

		// 等鸡啄完了米山、狗舔完了面山，火烧断了金锁
		group := sync.WaitGroup{}
		group.Add(5)
		go func() {
			for i := 0; i < len(invs); i++ {
				ctRe := <-fitChan.ctRe
				continuityList[ctRe.index] = ctRe.value
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				hpRe := <-fitChan.hpRe
				hairpinList[hpRe.index] = hpRe.value
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				for j := 0; j < len(seqSet); j++ {
					hmRe := <-fitChan.hmRe
					hmTable[hmRe.index1][hmRe.index2] = hmRe.value
				}
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				for j := 0; j < len(seqSet); j++ {
					if j == index {
						continue
					}
					smRe := <-fitChan.smRe
					smTable[smRe.index1][smRe.index2] = smRe.value
				}
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				mtRe := <-fitChan.mtRe
				mtList[mtRe.index] = mtRe.value
			}
			group.Done()
		}()
		group.Wait()

		hmList := sumTableRow(hmTable)
		smList := sumTableRow(smTable)
		// 连续与发卡使用实际值即可
		// 但是对于H-测度与相似性，一条DNA链的改变会影响其他链的值
		// 考虑使用
		// 1、七条的总值 || 平均值
		// 2、它对于其他DNA链在两个数值上的贡献度
		minCt := min(minCt, slices.Min(continuityList))
		minHp := min(minHp, slices.Min(hairpinList))
		minHm := min(minHm, slices.Min(hmList))
		minSm := min(minSm, slices.Min(smList))

		maxCt := slices.Max(continuityList)
		maxHp := slices.Max(hairpinList)
		maxHm := slices.Max(hmList)
		maxSm := slices.Max(smList)

		// Norm
		norm(continuityList, minCt, maxCt)
		norm(hairpinList, minHp, maxHp)
		norm(hmList, minHm, maxHm)
		norm(smList, minSm, maxSm)

		mtDiviant(mtList, mtValues)

		// Fitness
		//fitness := make([]float64, len(invs))
		bestInd := 0
		var bestFit float64 = 0
		for i := 0; i < len(invs); i++ {
			invs[i].ct = continuityList[i]
			invs[i].hp = hairpinList[i]
			invs[i].hm = hmList[i]
			invs[i].sm = smList[i]
			invs[i].mt = mtList[i]
			invs[i].fitness = continuityList[i] * hairpinList[i] * hmList[i] * smList[i] * mtList[i]
			if i == 0 || invs[i].fitness < bestFit {
				bestInd = i
				bestFit = invs[i].fitness
			}
		}

		return bestInd
	}
}

func mtDiviant(values, compared []float64) {
	var sumCompared float64 = 0
	for _, value := range compared {
		sumCompared += value
	}
	for i, value := range values {
		avg := (sumCompared + value) / float64(len(compared)+1)
		var div float64 = 0
		minVal := min(value, slices.Min(compared))
		maxVal := max(value, slices.Max(compared))
		for _, c := range compared {
			div += math.Pow((c-avg)/(maxVal-minVal), 2)
		}
		div += math.Pow((value-avg)/(maxVal-minVal), 2)
		div /= DNASIZE
		values[i] = max(div, MINVALUE)
	}
}

func sumTableRow[T int | float64](table [][]T) []float64 {
	sumList := make([]float64, len(table))
	for i, row := range table {
		var sum T = 0
		for _, element := range row {
			sum += element
		}
		sumList[i] = float64(sum)
	}
	return sumList
}

func avgTableRow[T int | float64](table [][]T) []float64 {
	avgs := make([]float64, len(table))
	for i, row := range table {
		var sum T = 0
		for _, element := range row {
			sum += element
		}
		avgs[i] = float64(sum) / float64(len(row))
	}
	return avgs
}

func norm(values []float64, minVal, maxVal float64) {
	for i, _ := range values {
		values[i] = (values[i] - minVal) / (maxVal - minVal)
		values[i] = max(values[i], MINVALUE)
	}
}

func sum(l []float64) float64 {
	var s float64
	for i := range l {
		s += l[i]
	}
	return s
}
