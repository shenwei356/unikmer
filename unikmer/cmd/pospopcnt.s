// +build gc

#include "textflag.h"

// https://stackoverflow.com/questions/63248047
// func Pospopcnt(counts *[8]int32, buf []byte)
TEXT ·Pospopcnt(SB),NOSPLIT,$0-32
	MOVQ counts+0(FP), DI
	MOVQ buf_base+8(FP), SI		// SI = &buf[0]
	MOVQ buf_len+16(FP), CX		// CX = len(buf)

	// load counts into register R8--R15
	MOVL 4*7(DI), R8
	MOVL 4*6(DI), R9
	MOVL 4*5(DI), R10
	MOVL 4*4(DI), R11
	MOVL 4*3(DI), R12
	MOVL 4*2(DI), R13
	MOVL 4*1(DI), R14
	MOVL 4*0(DI), R15

	SUBQ $32, CX			// pre-subtract 32 bit from CX
	JL scalar

vector:	VMOVDQU (SI), Y0		// load 32 bytes from buf
	ADDQ $32, SI			// advance SI past them

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R15			// add to counter
	VPADDD Y0, Y0, Y0		// shift Y0 left by one place

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R14			// add to counter
	VPADDD Y0, Y0, Y0		// shift Y0 left by one place

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R13			// add to counter
	VPADDD Y0, Y0, Y0		// shift Y0 left by one place

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R12			// add to counter
	VPADDD Y0, Y0, Y0		// shift Y0 left by one place

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R11			// add to counter
	VPADDD Y0, Y0, Y0		// shift Y0 left by one place

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R10			// add to counter
	VPADDD Y0, Y0, Y0		// shift Y0 left by one place

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R9			// add to counter
	VPADDD Y0, Y0, Y0		// shift Y0 left by one place

	VPMOVMSKB Y0, AX		// move MSB of Y0 bytes to AX
	POPCNTL AX, AX			// count population of AX
	ADDL AX, R8			// add to counter

	SUBQ $32, CX
	JGE vector			// repeat as long as bytes are left

scalar:	ADDQ $32, CX			// undo last subtraction
	JE done				// if CX=0, there's nothing left

loop:	MOVBLZX (SI), AX		// load a byte from buf
	INCQ SI				// advance past it

	BTL $0, AX			// is bit 0 set?
	ADCL $0, R8			// add it to R8

	BTL $1, AX			// is bit 1 set?
	ADCL $0, R9			// add it to R9

	BTL $2, AX			// is bit 2 set?
	ADCL $0, R10			// add it to R10

	BTL $3, AX			// is bit 3 set?
	ADCL $0, R11			// add it to R11

	BTL $4, AX			// is bit 4 set?
	ADCL $0, R12			// add it to R12

	BTL $5, AX			// is bit 5 set?
	ADCL $0, R13			// add it to R13

	BTL $6, AX			// is bit 6 set?
	ADCL $0, R14			// add it to R14

	BTL $7, AX			// is bit 7 set?
	ADCL $0, R15			// add it to R15

	DECQ CX				// mark this byte as done
	JNE loop			// and proceed if any bytes are left

	// write R8--R15 back to counts
done:	MOVL R8, 4*7(DI)
	MOVL R9, 4*6(DI)
	MOVL R10, 4*5(DI)
	MOVL R11, 4*4(DI)
	MOVL R12, 4*3(DI)
	MOVL R13, 4*2(DI)
	MOVL R14, 4*1(DI)
	MOVL R15, 4*0(DI)

	VZEROUPPER			// restore SSE-compatibility
	RET
