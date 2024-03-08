package GoDNA

import (
	"testing"
)

func TestToStr(t *testing.T) {
	seq := make(Seq, 6)
	seq = Seq([]int{0, 0, 1, 1, 2, 2, 3, 3, 6})
	str, err := seq.ToStr()
	if err == nil {
		t.Logf(str)
	} else {
		t.Log("error: ", err)
	}
	//fmt.Println(seq.ToStr())
}
