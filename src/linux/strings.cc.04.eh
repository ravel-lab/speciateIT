
;; Function char* strpre(const char*, const char*) (_Z6strprePKcS0_)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Forwarding edge 6->7 to 5 failed.
Merged 9 and 10 without moving.
Merged 9 and 11 without moving.


try_optimize_cfg iteration 2

Forwarding edge 6->7 to 5 failed.


try_optimize_cfg iteration 1

Forwarding edge 5->6 to 4 failed.
(note 1 0 12 ("./strings.cc") 43)

;; Start of basic block 0, registers live: (nil)
(note 12 1 8 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 8 12 9 0 ./strings.cc:43 (set (reg/v/f:DI 65 [ str ])
        (reg:DI 5 di [ str ])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./strings.cc:43 (set (reg/v/f:DI 66 [ pre ])
        (reg:DI 4 si [ pre ])) -1 (nil)
    (nil))

(note 10 9 15 0 NOTE_INSN_FUNCTION_BEG)

(note 15 10 16 0 ("./strings.cc") 47)

(insn 16 15 17 0 ./strings.cc:47 (set (reg:QI 60 [ temp.36 ])
        (mem:QI (reg/v/f:DI 66 [ pre ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 17 16 18 0 ./strings.cc:47 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ temp.36 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 18 17 20 0 ./strings.cc:47 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 29)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 20 18 21 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 21 20 22 1 ("./strings.cc") 48)

(insn 22 21 23 1 ./strings.cc:48 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ temp.36 ])
            (mem:QI (reg/v/f:DI 65 [ str ]) [0 S1 A8]))) -1 (nil)
    (nil))

(jump_insn 23 22 25 1 ./strings.cc:48 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 52)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 25 23 26 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 26 25 27 2 ./strings.cc:48 (set (reg:DI 61 [ ivtmp.33 ])
        (reg/v/f:DI 66 [ pre ])) -1 (nil)
    (nil))

(jump_insn 27 26 28 2 ./strings.cc:48 (set (pc)
        (label_ref 41)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

(barrier 28 27 29)

;; Start of basic block 3, registers live: (nil)
(code_label 29 28 30 3 3 "" [2 uses])

(note 30 29 31 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 31 30 32 3 ("./strings.cc") 47)

(insn 32 31 33 3 ./strings.cc:47 (set (reg:DI 62 [ D.2513 ])
        (reg/v/f:DI 65 [ str ])) -1 (nil)
    (nil))

(jump_insn 33 32 34 3 ./strings.cc:47 (set (pc)
        (label_ref 55)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 34 33 35)

;; Start of basic block 4, registers live: (nil)
(code_label 35 34 36 4 9 "" [1 uses])

(note 36 35 37 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 37 36 38 4 ./strings.cc:47 (parallel [
            (set (reg:DI 61 [ ivtmp.33 ])
                (plus:DI (reg:DI 61 [ ivtmp.33 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 38 37 39 4 ("./strings.cc") 48)

(insn 39 38 40 4 ./strings.cc:48 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 63 [ D.2512 ])
            (mem:QI (reg/v/f:DI 65 [ str ]) [0 S1 A8]))) -1 (nil)
    (nil))

(jump_insn 40 39 41 4 ./strings.cc:48 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 52)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 41 40 42 5 7 "" [1 uses])

(note 42 41 43 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 43 42 44 5 ("./strings.cc") 49)

(insn 44 43 45 5 ./strings.cc:49 (parallel [
            (set (reg/v/f:DI 65 [ str ])
                (plus:DI (reg/v/f:DI 65 [ str ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 45 44 46 5 ("./strings.cc") 47)

(insn 46 45 47 5 ./strings.cc:47 (set (reg:QI 63 [ D.2512 ])
        (mem:QI (plus:DI (reg:DI 61 [ ivtmp.33 ])
                (const_int 1 [0x1])) [0 S1 A8])) -1 (nil)
    (nil))

(insn 47 46 48 5 ./strings.cc:47 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 63 [ D.2512 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 48 47 51 5 ./strings.cc:47 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 29)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(note 51 48 49 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(jump_insn 49 51 50 6 ./strings.cc:47 (set (pc)
        (label_ref 35)) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

(barrier 50 49 52)

;; Start of basic block 7, registers live: (nil)
(code_label 52 50 53 7 5 "" [2 uses])

(note 53 52 54 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 54 53 55 7 ./strings.cc:47 (set (reg:DI 62 [ D.2513 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 55 54 56 8 8 "" [1 uses])

(note 56 55 57 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 57 56 60 8 ./strings.cc:47 (set (reg:DI 64 [ <result> ])
        (reg:DI 62 [ D.2513 ])) -1 (nil)
    (nil))

(note 60 57 61 8 NOTE_INSN_FUNCTION_END)

(note 61 60 63 8 ("./strings.cc") 51)

(insn 63 61 69 8 ./strings.cc:51 (set (reg/i:DI 0 ax)
        (reg:DI 64 [ <result> ])) -1 (nil)
    (nil))

(insn 69 63 0 8 ./strings.cc:51 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)


;; Function char* strlop(char*, char) (_Z6strlopPcc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Deleted label in block 5.
Forwarding edge 5->6 to 8 failed.
Merged 8 and 9 without moving.
Merged 8 and 10 without moving.


try_optimize_cfg iteration 2

Forwarding edge 5->6 to 8 failed.


try_optimize_cfg iteration 1

Forwarding edge 4->5 to 7 failed.
(note 1 0 8 ("./strings.cc") 81)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./strings.cc:81 (set (reg/v/f:DI 62 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./strings.cc:81 (set (reg:SI 64)
        (reg:SI 4 si [ ch ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./strings.cc:81 (set (reg/v:QI 63 [ ch ])
        (subreg:QI (reg:SI 64) 0)) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./strings.cc") 84)

(insn 11 10 12 0 ./strings.cc:84 (set (reg:QI 58 [ temp.67 ])
        (mem:QI (reg/v/f:DI 62 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./strings.cc:84 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 58 [ temp.67 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 13 12 15 0 ./strings.cc:84 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 36)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 15 13 16 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 16 15 17 1 ./strings.cc:84 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 58 [ temp.67 ])
            (reg/v:QI 63 [ ch ]))) -1 (nil)
    (nil))

(jump_insn 17 16 19 1 ./strings.cc:84 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 36)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 19 17 20 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 20 19 21 2 ./strings.cc:84 (set (reg/v/f:DI 59 [ last ])
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 21 20 22 3 19 "" [1 uses])

(note 22 21 23 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 23 22 24 3 ("./strings.cc") 85)

(insn 24 23 25 3 ./strings.cc:85 (parallel [
            (set (reg/v/f:DI 59 [ last ])
                (plus:DI (reg/v/f:DI 59 [ last ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 25 24 26 3 ("./strings.cc") 84)

(insn 26 25 27 3 ./strings.cc:84 (set (reg:QI 60 [ D.2553 ])
        (mem:QI (reg/v/f:DI 59 [ last ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 27 26 28 3 ./strings.cc:84 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ D.2553 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 28 27 30 3 ./strings.cc:84 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 39)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 30 28 31 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 31 30 32 4 ./strings.cc:84 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ D.2553 ])
            (reg/v:QI 63 [ ch ]))) -1 (nil)
    (nil))

(jump_insn 32 31 35 4 ./strings.cc:84 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 21)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 35 32 33 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 33 35 34 5 ./strings.cc:84 (set (pc)
        (label_ref 39)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 34 33 36)

;; Start of basic block 6, registers live: (nil)
(code_label 36 34 37 6 16 "" [2 uses])

(note 37 36 38 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 38 37 39 6 ./strings.cc:84 (set (reg/v/f:DI 59 [ last ])
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 39 38 40 7 20 "" [2 uses])

(note 40 39 41 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 41 40 42 7 ("./strings.cc") 86)

(insn 42 41 43 7 ./strings.cc:86 (set (mem:QI (reg/v/f:DI 59 [ last ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 43 42 46 7 ./strings.cc:86 (set (reg:DI 61 [ <result> ])
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))

(note 46 43 47 7 NOTE_INSN_FUNCTION_END)

(note 47 46 49 7 ("./strings.cc") 88)

(insn 49 47 55 7 ./strings.cc:88 (set (reg/i:DI 0 ax)
        (reg:DI 61 [ <result> ])) -1 (nil)
    (nil))

(insn 55 49 0 7 ./strings.cc:88 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)


;; Function char* strsplit(char*, char) (_Z8strsplitPcc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Forwarding edge 4->5 to 7 failed.
Deleted label in block 9.
Deleted label in block 10.
Deleted label in block 12.
Merged 13 and 14 without moving.
Merged 13 and 15 without moving.


try_optimize_cfg iteration 2

Forwarding edge 4->5 to 7 failed.


try_optimize_cfg iteration 1

Forwarding edge 3->4 to 6 failed.
(note 1 0 8 ("./strings.cc") 126)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./strings.cc:126 (set (reg/v/f:DI 63 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./strings.cc:126 (set (reg:SI 65)
        (reg:SI 4 si [ ch ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./strings.cc:126 (set (reg/v:QI 64 [ ch ])
        (subreg:QI (reg:SI 65) 0)) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./strings.cc") 128)

(insn 11 10 12 0 ./strings.cc:128 (set (reg:QI 61 [ D.2584 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./strings.cc:128 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 61 [ D.2584 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 13 12 15 0 ./strings.cc:128 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 19)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 6700 [0x1a2c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 15 13 16 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 16 15 17 1 ./strings.cc:128 (set (reg/v/f:DI 63 [ s ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 17 16 18 1 ./strings.cc:128 (set (pc)
        (label_ref 71)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 18 17 19)

;; Start of basic block 2, registers live: (nil)
(code_label 19 18 20 2 26 "" [1 uses])

(note 20 19 21 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 21 20 22 2 ("./strings.cc") 129)

(insn 22 21 23 2 ./strings.cc:129 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 61 [ D.2584 ])
            (reg/v:QI 64 [ ch ]))) -1 (nil)
    (nil))

(jump_insn 23 22 25 2 ./strings.cc:129 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 43)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 25 23 26 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 26 25 27 3 ./strings.cc:129 (parallel [
            (set (reg/v/f:DI 63 [ s ])
                (plus:DI (reg/v/f:DI 63 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 27 26 28 3 ./strings.cc:129 (set (reg:QI 59 [ temp.105 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 28 27 29 3 ./strings.cc:129 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 59 [ temp.105 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 29 28 32 3 ./strings.cc:129 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1100 [0x44c])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 32 29 30 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(jump_insn 30 32 31 4 ./strings.cc:129 (set (pc)
        (label_ref 39)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

(barrier 31 30 33)

;; Start of basic block 5, registers live: (nil)
(code_label 33 31 34 5 32 "" [1 uses])

(note 34 33 35 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 35 34 36 5 ./strings.cc:129 (parallel [
            (set (reg/v/f:DI 63 [ s ])
                (plus:DI (reg/v/f:DI 63 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 36 35 37 5 ./strings.cc:129 (set (reg:QI 59 [ temp.105 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 37 36 38 5 ./strings.cc:129 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 59 [ temp.105 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 38 37 39 5 ./strings.cc:129 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1100 [0x44c])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(code_label 39 38 40 6 31 "" [1 uses])

(note 40 39 41 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 41 40 42 6 ./strings.cc:129 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:QI 64 [ ch ])
            (reg:QI 59 [ temp.105 ]))) -1 (nil)
    (nil))

(jump_insn 42 41 43 6 ./strings.cc:129 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 33)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 43 42 44 7 29 "" [1 uses])

(note 44 43 45 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 45 44 46 7 ("./strings.cc") 130)

(insn 46 45 47 7 ./strings.cc:130 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:QI 64 [ ch ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 47 46 49 7 ./strings.cc:130 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3300 [0xce4])
        (nil)))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(note 49 47 50 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 50 49 51 8 ("./strings.cc") 132)

(insn 51 50 52 8 ./strings.cc:132 (set (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 52 51 53 8 ./strings.cc:132 (parallel [
            (set (reg/v/f:DI 63 [ s ])
                (plus:DI (reg/v/f:DI 63 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 53 52 54 8 ("./strings.cc") 133)

(insn 54 53 55 8 ./strings.cc:133 (set (reg:QI 58 [ temp.106 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 55 54 56 8 ./strings.cc:133 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 58 [ temp.106 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 56 55 58 8 ./strings.cc:133 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 58 56 59 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 59 58 60 9 ./strings.cc:133 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 58 [ temp.106 ])
            (reg/v:QI 64 [ ch ]))) -1 (nil)
    (nil))

(jump_insn 60 59 61 9 ./strings.cc:133 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(code_label 61 60 62 10 35 "" [1 uses])

(note 62 61 63 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(insn 63 62 64 10 ./strings.cc:133 (parallel [
            (set (reg/v/f:DI 63 [ s ])
                (plus:DI (reg/v/f:DI 63 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 64 63 65 10 ./strings.cc:133 (set (reg:QI 60 [ temp.104 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 65 64 66 10 ./strings.cc:133 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ temp.104 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 66 65 68 10 ./strings.cc:133 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 10, registers live:
 (nil)

;; Start of basic block 11, registers live: (nil)
(note 68 66 69 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 69 68 70 11 ./strings.cc:133 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ temp.104 ])
            (reg:QI 58 [ temp.106 ]))) -1 (nil)
    (nil))

(jump_insn 70 69 71 11 ./strings.cc:133 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 61)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(code_label 71 70 72 12 28 "" [7 uses])

(note 72 71 73 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(insn 73 72 76 12 ./strings.cc:133 (set (reg:DI 62 [ <result> ])
        (reg/v/f:DI 63 [ s ])) -1 (nil)
    (nil))

(note 76 73 77 12 NOTE_INSN_FUNCTION_END)

(note 77 76 79 12 ("./strings.cc") 136)

(insn 79 77 85 12 ./strings.cc:136 (set (reg/i:DI 0 ax)
        (reg:DI 62 [ <result> ])) -1 (nil)
    (nil))

(insn 85 79 0 12 ./strings.cc:136 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 12, registers live:
 (nil)


;; Function char* strtrimto(char*, char) (_Z9strtrimtoPcc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Merged 5 and 6 without moving.
Merged 5 and 7 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 8 ("./strings.cc") 155)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./strings.cc:155 (set (reg/v/f:DI 62 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./strings.cc:155 (set (reg:SI 64)
        (reg:SI 4 si [ ch ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./strings.cc:155 (set (reg/v:QI 63 [ ch ])
        (subreg:QI (reg:SI 64) 0)) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./strings.cc") 157)

(insn 11 10 12 0 ./strings.cc:157 (set (reg:SI 60 [ D.2619 ])
        (sign_extend:SI (reg/v:QI 63 [ ch ]))) -1 (nil)
    (nil))

(insn 12 11 13 0 ./strings.cc:157 (set (reg:SI 4 si)
        (reg:SI 60 [ D.2619 ])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./strings.cc:157 (set (reg:DI 5 di)
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))

(call_insn/u 14 13 15 0 ./strings.cc:157 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strrchr") [flags 0x41] <function_decl 0x2ad10a49e400 strrchr>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
                (nil)))))

(insn 15 14 16 0 ./strings.cc:157 (set (reg:DI 65)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 16 15 17 0 ./strings.cc:157 (set (reg/v/f:DI 59 [ last ])
        (reg:DI 65)) -1 (nil)
    (nil))

(note 17 16 18 0 ("./strings.cc") 158)

(insn 18 17 19 0 ./strings.cc:158 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 59 [ last ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 19 18 21 0 ./strings.cc:158 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 23)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8100 [0x1fa4])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 21 19 22 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 22 21 23 1 ./strings.cc:158 (set (reg/v/f:DI 59 [ last ])
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 23 22 24 2 45 "" [1 uses])

(note 24 23 25 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 25 24 26 2 ("./strings.cc") 159)

(insn 26 25 27 2 ./strings.cc:159 (set (mem:QI (reg/v/f:DI 59 [ last ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 27 26 28 2 ("./strings.cc") 161)

(insn 28 27 29 2 ./strings.cc:161 (set (reg:SI 4 si)
        (reg:SI 60 [ D.2619 ])) -1 (nil)
    (nil))

(insn 29 28 30 2 ./strings.cc:161 (set (reg:DI 5 di)
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))

(call_insn/u 30 29 31 2 ./strings.cc:161 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strchr") [flags 0x41] <function_decl 0x2ad10a498000 strchr>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
                (nil)))))

(insn 31 30 32 2 ./strings.cc:161 (set (reg:DI 66)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 32 31 33 2 ./strings.cc:161 (set (reg/v/f:DI 58 [ first ])
        (reg:DI 66)) -1 (nil)
    (nil))

(note 33 32 34 2 ("./strings.cc") 162)

(insn 34 33 35 2 ./strings.cc:162 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 58 [ first ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 35 34 37 2 ./strings.cc:162 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 39)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 6774 [0x1a76])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 37 35 38 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 38 37 39 3 ./strings.cc:162 (parallel [
            (set (reg/v/f:DI 58 [ first ])
                (plus:DI (reg/v/f:DI 59 [ last ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 39 38 40 4 47 "" [1 uses])

(note 40 39 41 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 41 40 42 4 ./strings.cc:162 (parallel [
            (set (reg:DI 67)
                (plus:DI (reg/v/f:DI 58 [ first ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 42 41 45 4 ./strings.cc:162 (set (reg:DI 61 [ <result> ])
        (reg:DI 67)) -1 (nil)
    (nil))

(note 45 42 46 4 NOTE_INSN_FUNCTION_END)

(note 46 45 48 4 ("./strings.cc") 164)

(insn 48 46 54 4 ./strings.cc:164 (set (reg/i:DI 0 ax)
        (reg:DI 61 [ <result> ])) -1 (nil)
    (nil))

(insn 54 48 0 4 ./strings.cc:164 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)


;; Function char* strchomp(char*) (_Z8strchompPc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Merged 8 and 9 without moving.
Merged 8 and 10 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("./strings.cc") 69)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./strings.cc:69 (set (reg/v/f:DI 63 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./strings.cc") 71)

(insn 9 8 10 0 ./strings.cc:71 (set (reg:DI 5 di)
        (reg/v/f:DI 63 [ s ])) -1 (nil)
    (nil))

(call_insn/u 10 9 11 0 ./strings.cc:71 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strlen") [flags 0x41] <function_decl 0x2ad10a49b200 strlen>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (nil))))

(insn 11 10 12 0 ./strings.cc:71 (set (reg:DI 65)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 12 11 13 0 ./strings.cc:71 (set (reg:DI 61 [ D.2533 ])
        (reg:DI 65)) -1 (nil)
    (nil))

(insn 13 12 14 0 ./strings.cc:71 (parallel [
            (set (reg/v/f:DI 60 [ last ])
                (plus:DI (reg:DI 61 [ D.2533 ])
                    (reg/v/f:DI 63 [ s ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 14 13 15 0 ("./strings.cc") 72)

(insn 15 14 16 0 ./strings.cc:72 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 63 [ s ])
            (reg/v/f:DI 60 [ last ]))) -1 (nil)
    (nil))

(jump_insn 16 15 18 0 ./strings.cc:72 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 45)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 18 16 19 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 19 18 20 1 ./strings.cc:72 (parallel [
            (set (reg/v/f:DI 58 [ last.167 ])
                (plus:DI (reg/v/f:DI 60 [ last ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 20 19 21 1 ./strings.cc:72 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem:QI (reg/v/f:DI 58 [ last.167 ]) [0 S1 A8])
            (const_int 10 [0xa]))) -1 (nil)
    (nil))

(jump_insn 21 20 23 1 ./strings.cc:72 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 27)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 23 21 24 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 24 23 25 2 ./strings.cc:72 (set (reg/v/f:DI 59 [ last.166 ])
        (reg/v/f:DI 58 [ last.167 ])) -1 (nil)
    (nil))

(jump_insn 25 24 26 2 ./strings.cc:72 (set (pc)
        (label_ref 48)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

(barrier 26 25 27)

;; Start of basic block 3, registers live: (nil)
(code_label 27 26 28 3 53 "" [1 uses])

(note 28 27 29 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 29 28 30 3 ./strings.cc:72 (set (reg/v/f:DI 59 [ last.166 ])
        (reg/v/f:DI 58 [ last.167 ])) -1 (nil)
    (nil))

(jump_insn 30 29 31 3 ./strings.cc:72 (set (pc)
        (label_ref 37)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 31 30 32)

;; Start of basic block 4, registers live: (nil)
(code_label 32 31 33 4 57 "" [1 uses])

(note 33 32 34 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 34 33 35 4 ./strings.cc:72 (parallel [
            (set (reg/v/f:DI 59 [ last.166 ])
                (plus:DI (reg/v/f:DI 59 [ last.166 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 35 34 36 4 ./strings.cc:72 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem:QI (reg/v/f:DI 59 [ last.166 ]) [0 S1 A8])
            (const_int 10 [0xa]))) -1 (nil)
    (nil))

(jump_insn 36 35 37 4 ./strings.cc:72 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 48)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 37 36 38 5 56 "" [1 uses])

(note 38 37 39 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 39 38 40 5 ./strings.cc:72 (parallel [
            (set (reg:DI 66)
                (minus:DI (reg/v/f:DI 60 [ last ])
                    (reg/v/f:DI 63 [ s ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 40 39 41 5 ./strings.cc:72 (set (reg:DI 68)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(insn 41 40 42 5 ./strings.cc:72 (parallel [
            (set (reg:DI 67)
                (minus:DI (reg:DI 68)
                    (reg:DI 66)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 42 41 43 5 ./strings.cc:72 (parallel [
            (set (reg:DI 69)
                (plus:DI (reg/v/f:DI 58 [ last.167 ])
                    (reg:DI 67)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 43 42 44 5 ./strings.cc:72 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 59 [ last.166 ])
            (reg:DI 69))) -1 (nil)
    (nil))

(jump_insn 44 43 45 5 ./strings.cc:72 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 32)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(code_label 45 44 46 6 51 "" [1 uses])

(note 46 45 47 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 47 46 48 6 ./strings.cc:72 (set (reg/v/f:DI 59 [ last.166 ])
        (reg/v/f:DI 63 [ s ])) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 48 47 49 7 55 "" [2 uses])

(note 49 48 50 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 50 49 51 7 ("./strings.cc") 74)

(insn 51 50 52 7 ./strings.cc:74 (set (mem:QI (plus:DI (reg/v/f:DI 59 [ last.166 ])
                (const_int 1 [0x1])) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 52 51 55 7 ./strings.cc:74 (set (reg:DI 62 [ <result> ])
        (reg/v/f:DI 63 [ s ])) -1 (nil)
    (nil))

(note 55 52 56 7 NOTE_INSN_FUNCTION_END)

(note 56 55 58 7 ("./strings.cc") 76)

(insn 58 56 64 7 ./strings.cc:76 (set (reg/i:DI 0 ax)
        (reg:DI 62 [ <result> ])) -1 (nil)
    (nil))

(insn 64 58 0 7 ./strings.cc:76 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)


;; Function char* strtrim(char*) (_Z7strtrimPc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Forwarding edge 0->2 to 4 failed.
Forwarding edge 0->2 to 4 failed.
Forwarding edge 5->6 to 8 failed.
Merged 9 and 10 without moving.
Merged 9 and 11 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->2 to 4 failed.
Forwarding edge 5->6 to 8 failed.


try_optimize_cfg iteration 1

Forwarding edge 0->1 to 3 failed.
Forwarding edge 4->5 to 7 failed.
(note 1 0 6 ("./strings.cc") 141)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./strings.cc:141 (set (reg/v/f:DI 64 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./strings.cc") 143)

(insn 9 8 10 0 ./strings.cc:143 (set (reg:DI 5 di)
        (reg/v/f:DI 64 [ s ])) -1 (nil)
    (nil))

(call_insn/u 10 9 11 0 ./strings.cc:143 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strlen") [flags 0x41] <function_decl 0x2ad10a49b200 strlen>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (nil))))

(insn 11 10 12 0 ./strings.cc:143 (set (reg:DI 66)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 12 11 13 0 ./strings.cc:143 (set (reg:DI 62 [ D.2596 ])
        (reg:DI 66)) -1 (nil)
    (nil))

(insn 13 12 14 0 ./strings.cc:143 (parallel [
            (set (reg:DI 67)
                (plus:DI (reg:DI 62 [ D.2596 ])
                    (reg/v/f:DI 64 [ s ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 14 13 15 0 ./strings.cc:143 (parallel [
            (set (reg/v/f:DI 58 [ last ])
                (plus:DI (reg:DI 67)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 15 14 16 0 ("./strings.cc") 144)

(insn 16 15 17 0 ./strings.cc:144 (set (reg:CC 17 flags)
        (compare:CC (reg/v/f:DI 64 [ s ])
            (reg/v/f:DI 58 [ last ]))) -1 (nil)
    (nil))

(jump_insn 17 16 20 0 ./strings.cc:144 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 36)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 20 17 18 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(jump_insn 18 20 19 1 ./strings.cc:144 (set (pc)
        (label_ref 28)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 19 18 21)

;; Start of basic block 2, registers live: (nil)
(code_label 21 19 22 2 90 "" [1 uses])

(note 22 21 23 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 23 22 24 2 ("./strings.cc") 145)

(insn 24 23 25 2 ./strings.cc:145 (parallel [
            (set (reg/v/f:DI 58 [ last ])
                (plus:DI (reg/v/f:DI 58 [ last ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 25 24 26 2 ("./strings.cc") 144)

(insn 26 25 27 2 ./strings.cc:144 (set (reg:CC 17 flags)
        (compare:CC (reg/v/f:DI 64 [ s ])
            (reg/v/f:DI 58 [ last ]))) -1 (nil)
    (nil))

(jump_insn 27 26 28 2 ./strings.cc:144 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 36)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 28 27 29 3 89 "" [1 uses])

(note 29 28 30 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 30 29 31 3 ./strings.cc:144 (set (reg:SI 68)
        (sign_extend:SI (mem:QI (reg/v/f:DI 58 [ last ]) [0 S1 A8]))) -1 (nil)
    (nil))

(insn 31 30 32 3 ./strings.cc:144 (set (reg:SI 5 di)
        (reg:SI 68)) -1 (nil)
    (nil))

(call_insn/u 32 31 33 3 ./strings.cc:144 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("isspace") [flags 0x41] <function_decl 0x2ad10a4b2900 isspace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 33 32 34 3 ./strings.cc:144 (set (reg:SI 61 [ D.2604 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 34 33 35 3 ./strings.cc:144 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 61 [ D.2604 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 35 34 36 3 ./strings.cc:144 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 21)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 36 35 37 4 88 "" [2 uses])

(note 37 36 38 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 38 37 39 4 ("./strings.cc") 146)

(insn 39 38 40 4 ./strings.cc:146 (set (mem:QI (plus:DI (reg/v/f:DI 58 [ last ])
                (const_int 1 [0x1])) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 40 39 41 4 ("./strings.cc") 147)

(insn 41 40 42 4 ./strings.cc:147 (set (reg:QI 60 [ D.2609 ])
        (mem:QI (reg/v/f:DI 64 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 42 41 43 4 ./strings.cc:147 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ D.2609 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 43 42 46 4 ./strings.cc:147 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 63)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 46 43 44 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 44 46 45 5 ./strings.cc:147 (set (pc)
        (label_ref 55)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 45 44 47)

;; Start of basic block 6, registers live: (nil)
(code_label 47 45 48 6 93 "" [1 uses])

(note 48 47 49 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 49 48 50 6 ("./strings.cc") 148)

(insn 50 49 51 6 ./strings.cc:148 (parallel [
            (set (reg/v/f:DI 64 [ s ])
                (plus:DI (reg/v/f:DI 64 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 51 50 52 6 ("./strings.cc") 147)

(insn 52 51 53 6 ./strings.cc:147 (set (reg:QI 60 [ D.2609 ])
        (mem:QI (reg/v/f:DI 64 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 53 52 54 6 ./strings.cc:147 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ D.2609 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 54 53 55 6 ./strings.cc:147 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 63)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 55 54 56 7 92 "" [1 uses])

(note 56 55 57 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 57 56 58 7 ./strings.cc:147 (set (reg:SI 69)
        (sign_extend:SI (reg:QI 60 [ D.2609 ]))) -1 (nil)
    (nil))

(insn 58 57 59 7 ./strings.cc:147 (set (reg:SI 5 di)
        (reg:SI 69)) -1 (nil)
    (nil))

(call_insn/u 59 58 60 7 ./strings.cc:147 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("isspace") [flags 0x41] <function_decl 0x2ad10a4b2900 isspace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 60 59 61 7 ./strings.cc:147 (set (reg:SI 59 [ D.2611 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 61 60 62 7 ./strings.cc:147 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 59 [ D.2611 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 62 61 63 7 ./strings.cc:147 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 47)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 63 62 64 8 91 "" [2 uses])

(note 64 63 65 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 65 64 68 8 ./strings.cc:147 (set (reg:DI 63 [ <result> ])
        (reg/v/f:DI 64 [ s ])) -1 (nil)
    (nil))

(note 68 65 69 8 NOTE_INSN_FUNCTION_END)

(note 69 68 71 8 ("./strings.cc") 150)

(insn 71 69 77 8 ./strings.cc:150 (set (reg/i:DI 0 ax)
        (reg:DI 63 [ <result> ])) -1 (nil)
    (nil))

(insn 77 71 0 8 ./strings.cc:150 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)


;; Function char* strsplit(char*) (_Z8strsplitPc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 5.
Forwarding edge 5->6 to 8 failed.
Merged 9 and 10 without moving.
Merged 9 and 11 without moving.


try_optimize_cfg iteration 2

Forwarding edge 5->6 to 8 failed.


try_optimize_cfg iteration 1

Forwarding edge 4->5 to 7 failed.
(note 1 0 6 ("./strings.cc") 108)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./strings.cc:108 (set (reg/v/f:DI 63 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./strings.cc") 110)

(insn 9 8 10 0 ./strings.cc:110 (set (reg:QI 61 [ D.2569 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 10 9 11 0 ./strings.cc:110 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 61 [ D.2569 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 ./strings.cc:110 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 24)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 6700 [0x1a2c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 13 11 14 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 14 13 15 1 ./strings.cc:110 (set (reg/v/f:DI 63 [ s ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 15 14 16 1 ./strings.cc:110 (set (pc)
        (label_ref 58)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 16 15 17)

;; Start of basic block 2, registers live: (nil)
(code_label 17 16 18 2 145 "" [1 uses])

(note 18 17 19 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 19 18 20 2 ("./strings.cc") 111)

(insn 20 19 21 2 ./strings.cc:111 (parallel [
            (set (reg/v/f:DI 63 [ s ])
                (plus:DI (reg/v/f:DI 63 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 21 20 22 2 ./strings.cc:111 (set (reg:QI 61 [ D.2569 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 22 21 23 2 ./strings.cc:111 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 61 [ D.2569 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 23 22 24 2 ./strings.cc:111 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 58)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1100 [0x44c])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 24 23 25 3 142 "" [1 uses])

(note 25 24 26 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 26 25 27 3 ./strings.cc:111 (set (reg:SI 64)
        (sign_extend:SI (reg:QI 61 [ D.2569 ]))) -1 (nil)
    (nil))

(insn 27 26 28 3 ./strings.cc:111 (set (reg:SI 5 di)
        (reg:SI 64)) -1 (nil)
    (nil))

(call_insn/u 28 27 29 3 ./strings.cc:111 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("isspace") [flags 0x41] <function_decl 0x2ad10a4b2900 isspace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 29 28 30 3 ./strings.cc:111 (set (reg:SI 60 [ D.2575 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 30 29 31 3 ./strings.cc:111 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 60 [ D.2575 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 31 30 33 3 ./strings.cc:111 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 17)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 33 31 34 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 34 33 35 4 ("./strings.cc") 114)

(insn 35 34 36 4 ./strings.cc:114 (set (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 36 35 37 4 ./strings.cc:114 (parallel [
            (set (reg/v/f:DI 63 [ s ])
                (plus:DI (reg/v/f:DI 63 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 37 36 38 4 ("./strings.cc") 115)

(insn 38 37 39 4 ./strings.cc:115 (set (reg:QI 58 [ temp.240 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 39 38 40 4 ./strings.cc:115 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 58 [ temp.240 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 40 39 43 4 ./strings.cc:115 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 58)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 43 40 41 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 41 43 42 5 ./strings.cc:115 (set (pc)
        (label_ref 50)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 42 41 44)

;; Start of basic block 6, registers live: (nil)
(code_label 44 42 45 6 148 "" [1 uses])

(note 45 44 46 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 46 45 47 6 ./strings.cc:115 (parallel [
            (set (reg/v/f:DI 63 [ s ])
                (plus:DI (reg/v/f:DI 63 [ s ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 47 46 48 6 ./strings.cc:115 (set (reg:QI 58 [ temp.240 ])
        (mem:QI (reg/v/f:DI 63 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 48 47 49 6 ./strings.cc:115 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 58 [ temp.240 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 49 48 50 6 ./strings.cc:115 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 58)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 50 49 51 7 147 "" [1 uses])

(note 51 50 52 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 52 51 53 7 ./strings.cc:115 (set (reg:SI 65)
        (sign_extend:SI (reg:QI 58 [ temp.240 ]))) -1 (nil)
    (nil))

(insn 53 52 54 7 ./strings.cc:115 (set (reg:SI 5 di)
        (reg:SI 65)) -1 (nil)
    (nil))

(call_insn/u 54 53 55 7 ./strings.cc:115 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("isspace") [flags 0x41] <function_decl 0x2ad10a4b2900 isspace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 55 54 56 7 ./strings.cc:115 (set (reg:SI 59 [ D.2579 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 56 55 57 7 ./strings.cc:115 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 59 [ D.2579 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 57 56 58 7 ./strings.cc:115 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 44)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 58 57 59 8 144 "" [4 uses])

(note 59 58 60 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 60 59 63 8 ./strings.cc:115 (set (reg:DI 62 [ <result> ])
        (reg/v/f:DI 63 [ s ])) -1 (nil)
    (nil))

(note 63 60 64 8 NOTE_INSN_FUNCTION_END)

(note 64 63 66 8 ("./strings.cc") 118)

(insn 66 64 72 8 ./strings.cc:118 (set (reg/i:DI 0 ax)
        (reg:DI 62 [ <result> ])) -1 (nil)
    (nil))

(insn 72 66 0 8 ./strings.cc:118 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)


;; Function char* strlopspace(char*) (_Z11strlopspacePc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 6 and 7 without moving.
Merged 6 and 8 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("./strings.cc") 93)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./strings.cc:93 (set (reg/v/f:DI 62 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./strings.cc") 96)

(insn 9 8 10 0 ./strings.cc:96 (set (reg:QI 60 [ D.2562 ])
        (mem:QI (reg/v/f:DI 62 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 10 9 11 0 ./strings.cc:96 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ D.2562 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 ./strings.cc:96 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 17)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 13 11 14 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 14 13 15 1 ./strings.cc:96 (set (reg/v/f:DI 58 [ last ])
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))

(jump_insn 15 14 16 1 ./strings.cc:96 (set (pc)
        (label_ref 38)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 16 15 17)

;; Start of basic block 2, registers live: (nil)
(code_label 17 16 18 2 155 "" [1 uses])

(note 18 17 19 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 19 18 20 2 ./strings.cc:96 (set (reg/v/f:DI 58 [ last ])
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))

(jump_insn 20 19 21 2 ./strings.cc:96 (set (pc)
        (label_ref 30)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

(barrier 21 20 22)

;; Start of basic block 3, registers live: (nil)
(code_label 22 21 23 3 159 "" [1 uses])

(note 23 22 24 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 24 23 25 3 ("./strings.cc") 97)

(insn 25 24 26 3 ./strings.cc:97 (parallel [
            (set (reg/v/f:DI 58 [ last ])
                (plus:DI (reg/v/f:DI 58 [ last ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 26 25 27 3 ("./strings.cc") 96)

(insn 27 26 28 3 ./strings.cc:96 (set (reg:QI 60 [ D.2562 ])
        (mem:QI (reg/v/f:DI 58 [ last ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 28 27 29 3 ./strings.cc:96 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 60 [ D.2562 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 29 28 30 3 ./strings.cc:96 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 38)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 30 29 31 4 158 "" [1 uses])

(note 31 30 32 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 32 31 33 4 ./strings.cc:96 (set (reg:SI 63)
        (sign_extend:SI (reg:QI 60 [ D.2562 ]))) -1 (nil)
    (nil))

(insn 33 32 34 4 ./strings.cc:96 (set (reg:SI 5 di)
        (reg:SI 63)) -1 (nil)
    (nil))

(call_insn/u 34 33 35 4 ./strings.cc:96 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("isspace") [flags 0x41] <function_decl 0x2ad10a4b2900 isspace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 35 34 36 4 ./strings.cc:96 (set (reg:SI 59 [ D.2564 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 36 35 37 4 ./strings.cc:96 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 59 [ D.2564 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 37 36 38 4 ./strings.cc:96 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 22)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 38 37 39 5 157 "" [2 uses])

(note 39 38 40 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 40 39 41 5 ("./strings.cc") 98)

(insn 41 40 42 5 ./strings.cc:98 (set (mem:QI (reg/v/f:DI 58 [ last ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 42 41 45 5 ./strings.cc:98 (set (reg:DI 61 [ <result> ])
        (reg/v/f:DI 62 [ s ])) -1 (nil)
    (nil))

(note 45 42 46 5 NOTE_INSN_FUNCTION_END)

(note 46 45 48 5 ("./strings.cc") 100)

(insn 48 46 54 5 ./strings.cc:100 (set (reg/i:DI 0 ax)
        (reg:DI 61 [ <result> ])) -1 (nil)
    (nil))

(insn 54 48 0 5 ./strings.cc:100 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)


;; Function char* strcasepre(const char*, const char*) (_Z10strcaseprePKcS0_)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 5.
Forwarding edge 5->6 to 4 failed.
Merged 8 and 9 without moving.
Merged 8 and 10 without moving.


try_optimize_cfg iteration 2

Forwarding edge 5->6 to 4 failed.


try_optimize_cfg iteration 1

Forwarding edge 4->5 to 3 failed.
(note 1 0 7 ("./strings.cc") 56)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./strings.cc:56 (set (reg/v/f:DI 64 [ str ])
        (reg:DI 5 di [ str ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./strings.cc:56 (set (reg/v/f:DI 65 [ pre ])
        (reg:DI 4 si [ pre ])) -1 (nil)
    (nil))

(note 5 4 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 5 11 0 ("./strings.cc") 60)

(insn 11 10 12 0 ./strings.cc:60 (set (reg:QI 62 [ D.2522 ])
        (mem:QI (reg/v/f:DI 65 [ pre ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./strings.cc:60 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 62 [ D.2522 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 13 12 14 0 ./strings.cc:60 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 19)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(code_label 14 13 15 1 164 "" [1 uses])

(note 15 14 16 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 16 15 17 1 ./strings.cc:60 (set (reg:DI 61 [ D.2523 ])
        (reg/v/f:DI 64 [ str ])) -1 (nil)
    (nil))

(jump_insn 17 16 18 1 ./strings.cc:60 (set (pc)
        (label_ref 50)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 18 17 19)

;; Start of basic block 2, registers live: (nil)
(code_label 19 18 20 2 163 "" [1 uses])

(note 20 19 21 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 21 20 22 2 ./strings.cc:60 (set (reg:DI 58 [ ivtmp.302 ])
        (reg/v/f:DI 65 [ pre ])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 22 21 23 3 166 "" [1 uses])

(note 23 22 24 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 24 23 25 3 ("./strings.cc") 61)

(insn 25 24 26 3 ./strings.cc:61 (set (reg:SI 66)
        (sign_extend:SI (reg:QI 62 [ D.2522 ]))) -1 (nil)
    (nil))

(insn 26 25 27 3 ./strings.cc:61 (set (reg:SI 5 di)
        (reg:SI 66)) -1 (nil)
    (nil))

(call_insn/u 27 26 28 3 ./strings.cc:61 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("tolower") [flags 0x41] <function_decl 0x2ad10a4b3500 tolower>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 28 27 29 3 ./strings.cc:61 (set (reg:SI 60 [ D.2525 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 29 28 30 3 ./strings.cc:61 (set (reg:SI 67)
        (sign_extend:SI (mem:QI (reg/v/f:DI 64 [ str ]) [0 S1 A8]))) -1 (nil)
    (nil))

(insn 30 29 31 3 ./strings.cc:61 (set (reg:SI 5 di)
        (reg:SI 67)) -1 (nil)
    (nil))

(call_insn/u 31 30 32 3 ./strings.cc:61 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("tolower") [flags 0x41] <function_decl 0x2ad10a4b3500 tolower>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 32 31 33 3 ./strings.cc:61 (set (reg:SI 59 [ D.2528 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 33 32 34 3 ./strings.cc:61 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 60 [ D.2525 ])
            (reg:SI 59 [ D.2528 ]))) -1 (nil)
    (nil))

(jump_insn 34 33 36 3 ./strings.cc:61 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 47)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 36 34 37 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 37 36 38 4 ("./strings.cc") 62)

(insn 38 37 39 4 ./strings.cc:62 (parallel [
            (set (reg/v/f:DI 64 [ str ])
                (plus:DI (reg/v/f:DI 64 [ str ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 39 38 40 4 ("./strings.cc") 60)

(insn 40 39 41 4 ./strings.cc:60 (set (reg:QI 62 [ D.2522 ])
        (mem:QI (plus:DI (reg:DI 58 [ ivtmp.302 ])
                (const_int 1 [0x1])) [0 S1 A8])) -1 (nil)
    (nil))

(insn 41 40 42 4 ./strings.cc:60 (parallel [
            (set (reg:DI 58 [ ivtmp.302 ])
                (plus:DI (reg:DI 58 [ ivtmp.302 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 42 41 43 4 ./strings.cc:60 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 62 [ D.2522 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 43 42 46 4 ./strings.cc:60 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 14)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 500 [0x1f4])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 46 43 44 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 44 46 45 5 ./strings.cc:60 (set (pc)
        (label_ref 22)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 45 44 47)

;; Start of basic block 6, registers live: (nil)
(code_label 47 45 48 6 167 "" [1 uses])

(note 48 47 49 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 49 48 50 6 ./strings.cc:60 (set (reg:DI 61 [ D.2523 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 50 49 51 7 165 "" [1 uses])

(note 51 50 52 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 52 51 55 7 ./strings.cc:60 (set (reg:DI 63 [ <result> ])
        (reg:DI 61 [ D.2523 ])) -1 (nil)
    (nil))

(note 55 52 56 7 NOTE_INSN_FUNCTION_END)

(note 56 55 58 7 ("./strings.cc") 64)

(insn 58 56 64 7 ./strings.cc:64 (set (reg/i:DI 0 ax)
        (reg:DI 63 [ <result> ])) -1 (nil)
    (nil))

(insn 64 58 0 7 ./strings.cc:64 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)


;; Function char* strclean(char*) (_Z8strcleanPc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 5.
Deleted label in block 6.
Merged 9 and 10 without moving.
Merged 9 and 11 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("./strings.cc") 29)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./strings.cc:29 (set (reg/v/f:DI 66 [ s ])
        (reg:DI 5 di [ s ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./strings.cc") 33)

(insn 9 8 10 0 ./strings.cc:33 (set (reg:QI 64 [ D.2495 ])
        (mem:QI (reg/v/f:DI 66 [ s ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 10 9 11 0 ./strings.cc:33 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 64 [ D.2495 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 ./strings.cc:33 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 17)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 13 11 14 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 14 13 15 1 ./strings.cc:33 (set (reg/v/f:DI 59 [ p ])
        (reg/v/f:DI 66 [ s ])) -1 (nil)
    (nil))

(jump_insn 15 14 16 1 ./strings.cc:33 (set (pc)
        (label_ref 56)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 16 15 17)

;; Start of basic block 2, registers live: (nil)
(code_label 17 16 18 2 172 "" [1 uses])

(note 18 17 19 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 19 18 20 2 ./strings.cc:33 (set (reg:DI 58 [ ivtmp.331 ])
        (reg/v/f:DI 66 [ s ])) -1 (nil)
    (nil))

(insn 20 19 21 2 ./strings.cc:33 (set (reg/v/f:DI 59 [ p ])
        (reg/v/f:DI 66 [ s ])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 21 20 22 3 175 "" [1 uses])

(note 22 21 23 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 23 22 24 3 ("./strings.cc") 34)

(insn 24 23 25 3 ./strings.cc:34 (set (reg:SI 63 [ D.2499 ])
        (sign_extend:SI (reg:QI 64 [ D.2495 ]))) -1 (nil)
    (nil))

(insn 25 24 26 3 ./strings.cc:34 (set (reg:SI 5 di)
        (reg:SI 63 [ D.2499 ])) -1 (nil)
    (nil))

(call_insn/u 26 25 27 3 ./strings.cc:34 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("isalnum") [flags 0x41] <function_decl 0x2ad10a4adb00 isalnum>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 27 26 28 3 ./strings.cc:34 (set (reg:SI 62 [ D.2500 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 28 27 29 3 ./strings.cc:34 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 62 [ D.2500 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 29 28 31 3 ./strings.cc:34 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 41)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3000 [0xbb8])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 31 29 32 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 32 31 33 4 ./strings.cc:34 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 64 [ D.2495 ])
            (const_int 45 [0x2d]))) -1 (nil)
    (nil))

(jump_insn 33 32 35 4 ./strings.cc:34 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 41)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 4880 [0x1310])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 35 33 36 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 36 35 37 5 ./strings.cc:34 (set (reg:SI 5 di)
        (reg:SI 63 [ D.2499 ])) -1 (nil)
    (nil))

(call_insn/u 37 36 38 5 ./strings.cc:34 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("isspace") [flags 0x41] <function_decl 0x2ad10a4b2900 isspace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 38 37 39 5 ./strings.cc:34 (set (reg:SI 61 [ D.2501 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 39 38 40 5 ./strings.cc:34 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 61 [ D.2501 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 40 39 41 5 ./strings.cc:34 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 49)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(code_label 41 40 42 6 176 "" [2 uses])

(note 42 41 43 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 43 42 44 6 ("./strings.cc") 35)

(insn 44 43 45 6 ./strings.cc:35 (set (reg:SI 5 di)
        (reg:SI 63 [ D.2499 ])) -1 (nil)
    (nil))

(call_insn/u 45 44 46 6 ./strings.cc:35 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("tolower") [flags 0x41] <function_decl 0x2ad10a4b3500 tolower>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (mem:BLK (scratch) [0 A8]))
        (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
            (nil))))

(insn 46 45 47 6 ./strings.cc:35 (set (reg:SI 60 [ D.2502 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 47 46 48 6 ./strings.cc:35 (set (mem:QI (reg/v/f:DI 59 [ p ]) [0 S1 A8])
        (subreg:QI (reg:SI 60 [ D.2502 ]) 0)) -1 (nil)
    (nil))

(insn 48 47 49 6 ./strings.cc:35 (parallel [
            (set (reg/v/f:DI 59 [ p ])
                (plus:DI (reg/v/f:DI 59 [ p ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 49 48 50 7 179 "" [1 uses])

(note 50 49 51 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 51 50 52 7 ("./strings.cc") 33)

(insn 52 51 53 7 ./strings.cc:33 (set (reg:QI 64 [ D.2495 ])
        (mem:QI (plus:DI (reg:DI 58 [ ivtmp.331 ])
                (const_int 1 [0x1])) [0 S1 A8])) -1 (nil)
    (nil))

(insn 53 52 54 7 ./strings.cc:33 (parallel [
            (set (reg:DI 58 [ ivtmp.331 ])
                (plus:DI (reg:DI 58 [ ivtmp.331 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 54 53 55 7 ./strings.cc:33 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 64 [ D.2495 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 55 54 56 7 ./strings.cc:33 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 21)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 56 55 57 8 174 "" [1 uses])

(note 57 56 58 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 58 57 59 8 ("./strings.cc") 36)

(insn 59 58 60 8 ./strings.cc:36 (set (mem:QI (reg/v/f:DI 59 [ p ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 60 59 63 8 ./strings.cc:36 (set (reg:DI 65 [ <result> ])
        (reg/v/f:DI 66 [ s ])) -1 (nil)
    (nil))

(note 63 60 64 8 NOTE_INSN_FUNCTION_END)

(note 64 63 66 8 ("./strings.cc") 38)

(insn 66 64 72 8 ./strings.cc:38 (set (reg/i:DI 0 ax)
        (reg:DI 65 [ <result> ])) -1 (nil)
    (nil))

(insn 72 66 0 8 ./strings.cc:38 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

