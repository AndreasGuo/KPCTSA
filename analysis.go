package GoDNA

import (
	"math"
	"slices"
)

func Continuity(seq Seq, t int) int {
	last := -1
	count := 1
	value := 0

	for _, base := range seq {
		if base == last {
			count += 1
		} else {
			if count >= t {
				value += count * count
			}
			count = 1
		}
		last = base
	}

	// continuous at tail
	if count >= t {
		value += count * count
	}

	return value
}

// define two bases are complementary or not
func cp(base1, base2 int) int {
	if base1+base2 == 3 {
		return 1
	}
	return 0
}

// threshold
func th(i, threshold int) int {
	if i > threshold {
		return i
	}
	return 0
}
func Hairpin(seq Seq, rMin int, pMin, pinLenT int) int {
	value := 0
	l := len(seq)
	// p is length of stem
	for p := pMin; p <= (l-rMin)/2; p++ {
		// r is length of loop
		for r := rMin; r <= l-2*p; r++ {
			// "i" is length of unused bases
			for i := 0; i <= l-2*p-r; i++ {
				sCp := 0
				// pin len defined shorter stem of a hairpin
				// suppose a real hair pin, it has a long side and a short side
				// which p+i is the long side
				pinLen := min(p+i, l-p-i-r)
				// check if complementary from both side of loop
				for j := 0; j < pinLen; j++ {
					sCp += cp(seq[p+i-j-1], seq[p+i+r+j])
				}
				value += th(sCp, pinLen/pinLenT)
			}
		}
	}
	return value
}

func shift(seq Seq, i int) Seq {
	// right shift i>0; left shift i<0
	l := len(seq)
	if i == 0 {
		return seq
	}
	absI := int(math.Abs(float64(i)))
	temp := make([]int, absI)
	for i := range temp {
		temp[i] = 5
	}
	switch {
	case 0 < i && i < l:
		return append(temp, seq[0:l-i]...)
	case i < 0 && i > -l:
		return append(seq[absI:l], temp...)
	default:
		temp = make([]int, l)
		for i := range temp {
			temp[i] = 5
		}
		return temp
	}
}

func lenNB(seq Seq) int {
	count := 0
	for _, base := range seq {
		if base <= G {
			count += 1
		}
	}
	return count
}

const HDIS = 0.17

func hDis(seq1, seq2 Seq) int {
	// this requires length of seq2 not less than seq1
	if len(seq2) < len(seq1) {
		panic("length of seq2 less than seq1 in hDis.")
	}
	cpS := 0
	l := len(seq1)
	for i := 0; i < l; i++ {
		cpS += cp(seq1[i], seq2[i])
	}
	threshold := HDIS * float64(lenNB(seq2)) / 2
	return th(cpS, int(threshold))
}

// continuous complementary start from i
func ccp(seq1, seq2 Seq, i int) int {
	c := 0
	for j := i; j < len(seq1); j++ {
		if seq1[j] == G-seq2[j] {
			c += 1
		} else {
			break
		}
	}
	return c
}

func hCon(seq1, seq2 Seq) int {
	HCON := 6
	value := 0
	for i := 0; i < len(seq1); i++ {
		value = value + th(ccp(seq1, seq2, i), HCON)
	}
	return value
}

func HMeasure(lSeq, rSeq Seq) int {
	l := len(lSeq)
	gap := int(math.Round(float64(l) / 4))
	HM := 0
	reversedLSeq := make(Seq, len(lSeq))
	copy(reversedLSeq, lSeq)
	slices.Reverse(reversedLSeq)
	for g := 0; g <= gap; g++ {
		temp := make([]int, g)
		for i := range temp {
			temp[i] = 5
		}
		tempSeq := append(rSeq, temp...)
		tempSeq = append(tempSeq, rSeq...)

		for i := -l + 1; i < l-1; i++ {
			shiftedSeq := shift(tempSeq, i)
			tempHM := hDis(reversedLSeq, shiftedSeq) + hCon(reversedLSeq, shiftedSeq)
			HM = max(HM, tempHM)
		}
	}
	return HM
}

const SDIS = 0.17

func sDis(seq1, seq2 Seq) int {
	sEq := 0
	for i := 0; i < len(seq1); i++ {
		if seq1[i] == seq2[i] {
			sEq += 1
		}
	}
	return th(sEq, int(SDIS*float64(len(seq2))/2))
}

func cEq(seq1, seq2 Seq, i int) int {
	c := 0
	for j := i; j < len(seq1); j++ {
		if seq1[j] == seq2[j] {
			c += 1
		} else {
			break
		}
	}
	return c
}

const SCON = 6

func sCon(seq1, seq2 Seq) int {
	value := 0
	for i := 0; i < len(seq1); i++ {
		value += th(cEq(seq1, seq2, i), SCON)
	}
	return value
}

func Similarity(lSeq, rSeq Seq) int {
	l := len(lSeq)
	gap := int(math.Round(float64(l) / 4))
	SM := 0
	//slices.Reverse(lSeq)
	for g := 0; g < gap; g++ {
		temp := make([]int, g)
		for i := range temp {
			temp[i] = 5
		}
		tempSeq := append(rSeq, temp...)
		tempSeq = append(tempSeq, rSeq...)

		for i := -l + 1; i < l; i++ {
			shiftedSeq := shift(tempSeq, i)
			tempSm := sDis(lSeq, shiftedSeq) + sCon(lSeq, shiftedSeq)
			SM = max(SM, tempSm)
		}
	}
	return SM
}
