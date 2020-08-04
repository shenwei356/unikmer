// +build gc

#include "textflag.h"

// https://stackoverflow.com/questions/63242918/golang-assembly-implement-of-mm-add-epi32
//
// func __mm_add_epi32(x, y [8]int32) [8]int32
TEXT ·__mm_add_epi32(SB),0,$0
    VMOVDQU x+0(FP), Y0
    VMOVDQU y+32(FP), Y1
    VPADDD  Y1, Y0, Y0
    VMOVDQU Y0, q+64(FP)
    VZEROUPPER
    RET
