
;; Function void __static_initialization_and_destruction_0(int, int) (_Z41__static_initialization_and_destruction_0ii)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Edge 0->4 redirected to 6
Deleted label in block 2.
Edge 2->4 redirected to 6
Deleted label in block 3.
Deleting block 4.
Deleting block 5.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 10 ("./DNAsequence.cc") 67)

;; Start of basic block 0, registers live: (nil)
(note 10 1 6 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 6 10 7 0 ./DNAsequence.cc:67 (set (reg/v:SI 58 [ __initialize_p ])
        (reg:SI 5 di [ __initialize_p ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./DNAsequence.cc:67 (set (reg/v:SI 59 [ __priority ])
        (reg:SI 4 si [ __priority ])) -1 (nil)
    (nil))

(note 8 7 12 0 NOTE_INSN_FUNCTION_BEG)

(note 12 8 13 0 ("./DNAsequence.cc") 67)

(insn 13 12 14 0 ./DNAsequence.cc:67 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 58 [ __initialize_p ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 14 13 16 0 ./DNAsequence.cc:67 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 36)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5467 [0x155b])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 16 14 17 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 17 16 18 1 ./DNAsequence.cc:67 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 59 [ __priority ])
            (const_int 65535 [0xffff]))) -1 (nil)
    (nil))

(jump_insn 18 17 20 1 ./DNAsequence.cc:67 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 36)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 7378 [0x1cd2])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 20 18 21 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 21 20 22 2 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

(insn 22 21 23 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2acec229c6e0 __ioinit>)) -1 (nil)
    (nil))

(call_insn 23 22 24 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitC1Ev") [flags 0x41] <function_decl 0x2acec1a5d400 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 24 23 25 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 1 dx)
        (symbol_ref:DI ("__dso_handle") [flags 0x40] <var_decl 0x2acec2433f20 __dso_handle>)) -1 (nil)
    (nil))

(insn 25 24 26 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 4 si)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 26 25 27 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("__tcf_0") [flags 0x3] <function_decl 0x2acec2427500 __tcf_0>)) -1 (nil)
    (nil))

(call_insn/j 27 26 28 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_atexit") [flags 0x41] <function_decl 0x2acec2427600 __cxa_atexit>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 2, registers live:
 (nil)

(barrier 28 27 31)

(note 31 28 32 NOTE_INSN_FUNCTION_END)

(note 32 31 36 ("./DNAsequence.cc") 67)

;; Start of basic block 3, registers live: (nil)
(code_label 36 32 39 3 5 "" [2 uses])

(note 39 36 0 3 [bb 3] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 3, registers live:
 (nil)


;; Function (static initializers for ./DNAsequence.cc) (_GLOBAL__I__Z11getAllKmersiRSt6vectorIPcSaIS0_EE)



try_optimize_cfg iteration 1

Deleting fallthru block 0.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 3 ("./DNAsequence.cc") 68)

(note 3 1 6 0 NOTE_INSN_FUNCTION_BEG)

;; Start of basic block 0, registers live: (nil)
(note 6 3 7 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(note 7 6 8 0 ("./DNAsequence.cc") 68)

(insn 8 7 9 0 ./DNAsequence.cc:68 (set (reg:SI 4 si)
        (const_int 65535 [0xffff])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./DNAsequence.cc:68 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 ./DNAsequence.cc:68 (call (mem:QI (symbol_ref:DI ("_Z41__static_initialization_and_destruction_0ii") [flags 0x3] <function_decl 0x2acec2427400 __static_initialization_and_destruction_0>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("./DNAsequence.cc") 68)


;; Function void __tcf_0(void*) (__tcf_0)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg/v/f:DI 58 [ D.36018 ])
        (reg:DI 5 di [ D.36018 ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

(insn 9 8 10 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2acec229c6e0 __ioinit>)) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitD1Ev") [flags 0x41] <function_decl 0x2acec1a5d600 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)


;; Function void std::vector<_Tp, _Alloc>::_M_insert_aux(__gnu_cxx::__normal_iterator<typename std::_Vector_base<_Tp, _Alloc>::_Tp_alloc_type::pointer, std::vector<_Tp, _Alloc> >, const _Tp&) [with _Tp = char*, _Alloc = std::allocator<char*>] (_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Redirecting jump 46 from 20 to 22.
Deleted label in block 6.
Deleted label in block 8.
Deleted label in block 11.
Deleted label in block 14.
Deleted label in block 16.
Deleted label in block 18.
Merged 19 and 20 without moving.
Merged 19 and 21 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 8 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 249)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:249 (set (reg/f:DI 80 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(insn 4 3 5 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:249 (set (reg/v:DI 81 [ __position ])
        (reg:DI 4 si [ __position ])) -1 (nil)
    (nil))

(insn 5 4 6 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:249 (set (reg/v/f:DI 82 [ __x ])
        (reg:DI 1 dx [ __x ])) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(insn 10 6 11 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:249 (set (reg:DI 61 [ __position$_M_current ])
        (reg/v:DI 81 [ __position ])) -1 (nil)
    (nil))

(note 11 10 12 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 251)

(insn 12 11 13 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:251 (set (reg:DI 79 [ D.35436 ])
        (mem/s/f:DI (plus:DI (reg/f:DI 80 [ this ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 13 12 14 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:251 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 79 [ D.35436 ])
            (mem/s/f:DI (plus:DI (reg/f:DI 80 [ this ])
                    (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 14 13 16 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:251 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 48)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 16 14 17 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 17 16 18 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:251 (set (reg/v/f:DI 73 [ __p ])
        (reg:DI 79 [ D.35436 ])) -1 (nil)
    (nil))

(note 18 17 19 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 19 18 20 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 73 [ __p ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 20 19 22 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 25)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 22 20 23 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 23 22 24 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:DI 83)
        (mem/f:DI (plus:DI (reg:DI 79 [ D.35436 ])
                (const_int -8 [0xfffffffffffffff8])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 24 23 25 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 73 [ __p ]) [5 S8 A64])
        (reg:DI 83)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 25 24 26 3 15 "" [1 uses])

(note 26 25 27 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 27 26 28 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 255)

(insn 28 27 29 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:255 (parallel [
            (set (reg:DI 58 [ temp.447 ])
                (plus:DI (reg:DI 79 [ D.35436 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 29 28 30 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:255 (set (mem/s:DI (plus:DI (reg/f:DI 80 [ this ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
        (reg:DI 58 [ temp.447 ])) -1 (nil)
    (nil))

(note 30 29 31 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 256)

(insn 31 30 32 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:256 (set (reg/v/f:DI 76 [ __x_copy ])
        (mem/f:DI (reg/v/f:DI 82 [ __x ]) [5 S8 A64])) -1 (nil)
    (nil))

(insn 32 31 33 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:256 (set (reg/v/f:DI 70 [ __first ])
        (reg:DI 61 [ __position$_M_current ])) -1 (nil)
    (nil))

(note 33 32 34 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h") 425)

(insn 34 33 35 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (parallel [
            (set (reg:DI 84)
                (plus:DI (reg:DI 58 [ temp.447 ])
                    (const_int -16 [0xfffffffffffffff0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 35 34 36 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (parallel [
            (set (reg:DI 85)
                (minus:DI (reg:DI 84)
                    (reg/v/f:DI 70 [ __first ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 36 35 37 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (parallel [
            (set (reg:DI 87)
                (ashiftrt:DI (reg:DI 85)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 85)
            (const_int 8 [0x8]))
        (nil)))

(insn 37 36 38 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (parallel [
            (set (reg:DI 71 [ D.37593 ])
                (ashift:DI (reg:DI 87)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 38 37 39 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (parallel [
            (set (reg:DI 88)
                (plus:DI (reg:DI 58 [ temp.447 ])
                    (const_int -8 [0xfffffffffffffff8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 39 38 40 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (parallel [
            (set (reg:DI 89)
                (minus:DI (reg:DI 88)
                    (reg:DI 71 [ D.37593 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 40 39 41 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (set (reg:DI 1 dx)
        (reg:DI 71 [ D.37593 ])) -1 (nil)
    (nil))

(insn 41 40 42 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (set (reg:DI 4 si)
        (reg/v/f:DI 70 [ __first ])) -1 (nil)
    (nil))

(insn 42 41 43 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (set (reg:DI 5 di)
        (reg:DI 89)) -1 (nil)
    (nil))

(call_insn 43 42 44 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:425 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memmove") [flags 0x41] <function_decl 0x2acebfb0b800 memmove>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(note 44 43 45 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 260)

(insn 45 44 46 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:260 (set (mem/f:DI (reg:DI 61 [ __position$_M_current ]) [5 S8 A64])
        (reg/v/f:DI 76 [ __x_copy ])) -1 (nil)
    (nil))

(jump_insn 46 45 47 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:260 (set (pc)
        (label_ref:DI 170)) 548 {jump} (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 47 46 48)

;; Start of basic block 4, registers live: (nil)
(code_label 48 47 49 4 13 "" [1 uses])

(note 49 48 50 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 50 49 51 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 809)

(insn 51 50 52 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (parallel [
            (set (reg:DI 90)
                (minus:DI (reg:DI 79 [ D.35436 ])
                    (mem/s/f:DI (reg/f:DI 80 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 52 51 53 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (parallel [
            (set (reg:DI 91)
                (ashiftrt:DI (reg:DI 90)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 90)
            (const_int 8 [0x8]))
        (nil)))

(insn 53 52 54 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 69 [ D.37630 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(note 54 53 55 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 265)

(insn 55 54 56 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:265 (set (reg:DI 92)
        (const_int 2305843009213693951 [0x1fffffffffffffff])) -1 (nil)
    (nil))

(insn 56 55 57 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:265 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 69 [ D.37630 ])
            (reg:DI 92))) -1 (nil)
    (nil))

(jump_insn 57 56 59 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:265 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 64)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 59 57 60 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 60 59 61 5 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 266)

(insn 61 60 62 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:266 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC0") [flags 0x2] <string_cst 0x2acec2675680>)) -1 (nil)
    (nil))

(call_insn 62 61 63 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:266 (call (mem:QI (symbol_ref:DI ("_ZSt20__throw_length_errorPKc") [flags 0x41] <function_decl 0x2acebfeb6a00 __throw_length_error>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 5, registers live:
 (nil)

(barrier 63 62 64)

;; Start of basic block 6, registers live: (nil)
(code_label 64 63 65 6 18 "" [1 uses])

(note 65 64 66 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 66 65 67 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 402)

(insn 67 66 68 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (set (reg/v:DI 75 [ __old_size ])
        (reg:DI 69 [ D.37630 ])) -1 (nil)
    (nil))

(note 68 67 69 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 271)

(insn 69 68 70 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:271 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 69 [ D.37630 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 70 69 72 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:271 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 76)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 72 70 73 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 73 72 74 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:271 (set (reg/v:DI 74 [ __len ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(jump_insn 74 73 75 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:271 (set (pc)
        (label_ref 79)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

(barrier 75 74 76)

;; Start of basic block 8, registers live: (nil)
(code_label 76 75 77 8 20 "" [1 uses])

(note 77 76 78 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 78 77 79 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:271 (parallel [
            (set (reg/v:DI 74 [ __len ])
                (ashift:DI (reg/v:DI 75 [ __old_size ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(code_label 79 78 80 9 22 "" [1 uses])

(note 80 79 81 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(note 81 80 82 9 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 272)

(insn 82 81 83 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 74 [ __len ])
            (reg/v:DI 75 [ __old_size ]))) -1 (nil)
    (nil))

(jump_insn 83 82 85 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (pc)
        (if_then_else (geu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 90)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(note 85 83 86 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(insn 86 85 87 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (reg/v:DI 74 [ __len ])
        (const_int 2305843009213693951 [0x1fffffffffffffff])) -1 (nil)
    (nil))

(insn 87 86 88 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (reg:DI 59 [ prephitmp.446 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 88 87 89 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (pc)
        (label_ref 96)) -1 (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 89 88 90)

;; Start of basic block 11, registers live: (nil)
(code_label 90 89 91 11 23 "" [1 uses])

(note 91 90 92 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 92 91 93 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (reg:DI 93)
        (const_int 2305843009213693951 [0x1fffffffffffffff])) -1 (nil)
    (nil))

(insn 93 92 94 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 74 [ __len ])
            (reg:DI 93))) -1 (nil)
    (nil))

(insn 94 93 95 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (reg:QI 94)
        (gtu:QI (reg:CC 17 flags)
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EQUAL (gtu:QI (reg/v:DI 74 [ __len ])
            (reg:DI 93))
        (nil)))

(insn 95 94 96 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:272 (set (reg:DI 59 [ prephitmp.446 ])
        (zero_extend:DI (reg:QI 94))) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(code_label 96 95 97 12 25 "" [1 uses])

(note 97 96 98 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(note 98 97 99 12 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 85)

(insn 99 98 100 12 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:85 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 59 [ prephitmp.446 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 100 99 102 12 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:85 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 106)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 12, registers live:
 (nil)

;; Start of basic block 13, registers live: (nil)
(note 102 100 103 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(note 103 102 104 13 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 86)

(call_insn 104 103 105 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:86 (call (mem:QI (symbol_ref:DI ("_ZSt17__throw_bad_allocv") [flags 0x41] <function_decl 0x2acebfeb6400 __throw_bad_alloc>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (nil))
;; End of basic block 13, registers live:
 (nil)

(barrier 105 104 106)

;; Start of basic block 14, registers live: (nil)
(code_label 106 105 107 14 26 "" [1 uses])

(note 107 106 108 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(note 108 107 109 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 88)

(insn 109 108 110 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (parallel [
            (set (reg:DI 67 [ D.37649 ])
                (ashift:DI (reg/v:DI 74 [ __len ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 110 109 111 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg:DI 5 di)
        (reg:DI 67 [ D.37649 ])) -1 (nil)
    (nil))

(call_insn 111 110 112 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Znwm") [flags 0x41] <function_decl 0x2acebfb9a700 operator new>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 112 111 113 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg:DI 68 [ D.37650 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 113 112 114 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg:DI 78 [ D.35410 ])
        (reg:DI 68 [ D.37650 ])) -1 (nil)
    (nil))

(insn 114 113 115 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg/v/f:DI 63 [ __first ])
        (mem/s/f:DI (reg/f:DI 80 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 115 114 116 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h") 300)

(insn 116 115 117 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 64 [ __last.112 ])
        (reg:DI 61 [ __position$_M_current ])) -1 (nil)
    (nil))

(insn 117 116 118 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 65 [ D.37714 ])
                (minus:DI (reg:DI 64 [ __last.112 ])
                    (reg/v/f:DI 63 [ __first ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 118 117 119 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 1 dx)
        (reg:DI 65 [ D.37714 ])) -1 (nil)
    (nil))

(insn 119 118 120 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 4 si)
        (reg/v/f:DI 63 [ __first ])) -1 (nil)
    (nil))

(insn 120 119 121 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 5 di)
        (reg:DI 78 [ D.35410 ])) -1 (nil)
    (nil))

(call_insn 121 120 122 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memmove") [flags 0x41] <function_decl 0x2acebfb0b800 memmove>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(note 122 121 123 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h") 301)

(insn 123 122 124 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:301 (parallel [
            (set (reg:DI 66 [ D.37701 ])
                (plus:DI (reg:DI 78 [ D.35410 ])
                    (reg:DI 65 [ D.37714 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 124 123 125 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:301 (set (reg/v/f:DI 72 [ __p ])
        (reg:DI 66 [ D.37701 ])) -1 (nil)
    (nil))

(note 125 124 126 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 126 125 127 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 72 [ __p ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 127 126 129 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 132)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 14, registers live:
 (nil)

;; Start of basic block 15, registers live: (nil)
(note 129 127 130 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(insn 130 129 131 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:DI 95)
        (mem/f:DI (reg/v/f:DI 82 [ __x ]) [5 S8 A64])) -1 (nil)
    (nil))

(insn 131 130 132 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 72 [ __p ]) [5 S8 A64])
        (reg:DI 95)) -1 (nil)
    (nil))
;; End of basic block 15, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(code_label 132 131 133 16 28 "" [1 uses])

(note 133 132 134 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(note 134 133 135 16 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 673)

(insn 135 134 136 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:673 (parallel [
            (set (reg:DI 60 [ __new_finish$_M_current ])
                (plus:DI (reg:DI 66 [ D.37701 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 136 135 137 16 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h") 300)

(insn 137 136 138 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 96)
        (mem/s/f:DI (plus:DI (reg/f:DI 80 [ this ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 138 137 139 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 62 [ D.37779 ])
                (minus:DI (reg:DI 96)
                    (reg:DI 64 [ __last.112 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (minus:DI (mem/s/f:DI (plus:DI (reg/f:DI 80 [ this ])
                    (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
            (reg:DI 64 [ __last.112 ]))
        (nil)))

(insn 139 138 140 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 1 dx)
        (reg:DI 62 [ D.37779 ])) -1 (nil)
    (nil))

(insn 140 139 141 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 4 si)
        (reg:DI 61 [ __position$_M_current ])) -1 (nil)
    (nil))

(insn 141 140 142 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 5 di)
        (reg:DI 60 [ __new_finish$_M_current ])) -1 (nil)
    (nil))

(call_insn 142 141 143 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memmove") [flags 0x41] <function_decl 0x2acebfb0b800 memmove>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(note 143 142 144 16 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 299)

(insn 144 143 145 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:299 (set (reg:DI 77 [ D.35467 ])
        (mem/s/f:DI (reg/f:DI 80 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 145 144 146 16 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 132)

(insn 146 145 147 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 77 [ D.35467 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 147 146 149 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 153)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3537 [0xdd1])
        (nil)))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(note 149 147 150 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(note 150 149 151 17 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 94)

(insn 151 150 152 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg:DI 77 [ D.35467 ])) -1 (nil)
    (nil))

(call_insn 152 151 153 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_ZdlPv") [flags 0x41] <function_decl 0x2acebfb9a900 operator delete>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(code_label 153 152 154 18 30 "" [1 uses])

(note 154 153 155 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(note 155 154 156 18 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 302)

(insn 156 155 157 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:302 (set (mem/s:DI (reg/f:DI 80 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])
        (reg:DI 78 [ D.35410 ])) -1 (nil)
    (nil))

(note 157 156 158 18 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 303)

(insn 158 157 159 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:303 (parallel [
            (set (reg:DI 97)
                (plus:DI (reg:DI 60 [ __new_finish$_M_current ])
                    (reg:DI 62 [ D.37779 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 159 158 160 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:303 (set (mem/s:DI (plus:DI (reg/f:DI 80 [ this ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
        (reg:DI 97)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (reg:DI 60 [ __new_finish$_M_current ])
            (reg:DI 62 [ D.37779 ]))
        (nil)))

(note 160 159 161 18 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 304)

(insn 161 160 162 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:304 (parallel [
            (set (reg:DI 98)
                (plus:DI (reg:DI 78 [ D.35410 ])
                    (reg:DI 67 [ D.37649 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 162 161 165 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:304 (set (mem/s:DI (plus:DI (reg/f:DI 80 [ this ])
                (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64])
        (reg:DI 98)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (reg:DI 78 [ D.35410 ])
            (reg:DI 67 [ D.37649 ]))
        (nil)))
;; End of basic block 18, registers live:
 (nil)

(note 165 162 166 NOTE_INSN_FUNCTION_END)

(note 166 165 170 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 304)

;; Start of basic block 19, registers live: (nil)
(code_label 170 166 173 19 32 "" [1 uses])

(note 173 170 0 19 [bb 19] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 19, registers live:
 (nil)


;; Function std::vector<_Tp, _Alloc>& std::vector<_Tp, _Alloc>::operator=(const std::vector<_Tp, _Alloc>&) [with _Tp = char*, _Alloc = std::allocator<char*>] (_ZNSt6vectorIPcSaIS0_EEaSERKS2_)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Deleted label in block 4.
Deleted label in block 6.
Deleted label in block 9.
Merged 12 and 13 without moving.
Merged 12 and 14 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 133)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:133 (set (reg/f:DI 74 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(insn 4 3 5 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:133 (set (reg/v/f:DI 75 [ __x ])
        (reg:DI 4 si [ __x ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 135)

(insn 10 9 11 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:135 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 75 [ __x ])
            (reg/f:DI 74 [ this ]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:135 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 127)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1036 [0x40c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 13 11 14 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 14 13 15 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 342)

(insn 15 14 16 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:342 (set (reg:DI 65 [ D.38311 ])
        (mem/s/f:DI (reg/v/f:DI 75 [ __x ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 16 15 17 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 809)

(insn 17 16 18 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 67 [ D.38332 ])
        (reg:DI 65 [ D.38311 ])) -1 (nil)
    (nil))

(insn 18 17 19 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 76)
        (mem/s/f:DI (plus:DI (reg/v/f:DI 75 [ __x ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 19 18 20 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (parallel [
            (set (reg:DI 66 [ D.38333 ])
                (minus:DI (reg:DI 76)
                    (reg:DI 67 [ D.38332 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (minus:DI (mem/s/f:DI (plus:DI (reg/v/f:DI 75 [ __x ])
                    (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
            (reg:DI 67 [ D.38332 ]))
        (nil)))

(note 20 19 21 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 402)

(insn 21 20 22 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (parallel [
            (set (reg:DI 77)
                (ashiftrt:DI (reg:DI 66 [ D.38333 ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 66 [ D.38333 ])
            (const_int 8 [0x8]))
        (nil)))

(insn 22 21 23 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (set (reg/v:DI 71 [ __xlen ])
        (reg:DI 77)) -1 (nil)
    (nil))

(insn 23 22 24 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (set (reg/f:DI 63 [ this ])
        (reg/f:DI 74 [ this ])) -1 (nil)
    (nil))

(note 24 23 25 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 342)

(insn 25 24 26 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:342 (set (reg:DI 62 [ D.38345 ])
        (mem/s/f:DI (reg/f:DI 63 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 26 25 27 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 809)

(insn 27 26 28 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 64 [ D.38363 ])
        (reg:DI 62 [ D.38345 ])) -1 (nil)
    (nil))

(note 28 27 29 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 138)

(insn 29 28 30 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:138 (set (reg:DI 79)
        (mem/s/f:DI (plus:DI (reg/f:DI 63 [ this ])
                (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64])) -1 (nil)
    (nil))

(insn 30 29 31 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:138 (parallel [
            (set (reg:DI 78)
                (minus:DI (reg:DI 79)
                    (reg:DI 64 [ D.38363 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (minus:DI (mem/s/f:DI (plus:DI (reg/f:DI 63 [ this ])
                    (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64])
            (reg:DI 64 [ D.38363 ]))
        (nil)))

(insn 31 30 32 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:138 (parallel [
            (set (reg:DI 81)
                (ashiftrt:DI (reg:DI 78)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 78)
            (const_int 8 [0x8]))
        (nil)))

(insn 32 31 33 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:138 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 71 [ __xlen ])
            (reg:DI 81))) -1 (nil)
    (nil))

(jump_insn 33 32 35 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:138 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 78)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7000 [0x1b58])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 35 33 36 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 36 35 37 2 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 85)

(insn 37 36 38 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:85 (set (reg:DI 82)
        (const_int 2305843009213693951 [0x1fffffffffffffff])) -1 (nil)
    (nil))

(insn 38 37 39 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:85 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 71 [ __xlen ])
            (reg:DI 82))) -1 (nil)
    (nil))

(jump_insn 39 38 41 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:85 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 45)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 41 39 42 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 42 41 43 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 86)

(call_insn 43 42 44 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:86 (call (mem:QI (symbol_ref:DI ("_ZSt17__throw_bad_allocv") [flags 0x41] <function_decl 0x2acebfeb6400 __throw_bad_alloc>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 44 43 45)

;; Start of basic block 4, registers live: (nil)
(code_label 45 44 46 4 41 "" [1 uses])

(note 46 45 47 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 47 46 48 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 88)

(insn 48 47 49 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (parallel [
            (set (reg:DI 60 [ D.38432 ])
                (ashift:DI (reg/v:DI 71 [ __xlen ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 49 48 50 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg:DI 5 di)
        (reg:DI 60 [ D.38432 ])) -1 (nil)
    (nil))

(call_insn 50 49 51 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Znwm") [flags 0x41] <function_decl 0x2acebfb9a700 operator new>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 51 50 52 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg:DI 61 [ D.38433 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 52 51 53 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:88 (set (reg/v/f:DI 70 [ __tmp ])
        (reg:DI 61 [ D.38433 ])) -1 (nil)
    (nil))

(note 53 52 54 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h") 300)

(insn 54 53 55 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 1 dx)
        (reg:DI 66 [ D.38333 ])) -1 (nil)
    (nil))

(insn 55 54 56 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 4 si)
        (reg:DI 65 [ D.38311 ])) -1 (nil)
    (nil))

(insn 56 55 57 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 5 di)
        (reg/v/f:DI 70 [ __tmp ])) -1 (nil)
    (nil))

(call_insn 57 56 58 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memmove") [flags 0x41] <function_decl 0x2acebfb0b800 memmove>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(note 58 57 59 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 142)

(insn 59 58 60 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:142 (set (reg/v/f:DI 69 [ D.36577 ])
        (mem/s/f:DI (reg/f:DI 74 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 60 59 61 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 132)

(insn 61 60 62 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 69 [ D.36577 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 62 61 64 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 68)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3537 [0xdd1])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 64 62 65 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 65 64 66 5 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 94)

(insn 66 65 67 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg/v/f:DI 69 [ D.36577 ])) -1 (nil)
    (nil))

(call_insn 67 66 68 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_ZdlPv") [flags 0x41] <function_decl 0x2acebfb9a900 operator delete>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(code_label 68 67 69 6 43 "" [1 uses])

(note 69 68 70 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 70 69 71 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 147)

(insn 71 70 72 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:147 (set (mem/s:DI (reg/f:DI 74 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])
        (reg/v/f:DI 70 [ __tmp ])) -1 (nil)
    (nil))

(note 72 71 73 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 148)

(insn 73 72 74 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:148 (set (reg:DI 58 [ prephitmp.771 ])
        (reg:DI 60 [ D.38432 ])) -1 (nil)
    (nil))

(insn 74 73 75 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:148 (parallel [
            (set (reg:DI 83)
                (plus:DI (reg/v/f:DI 70 [ __tmp ])
                    (reg:DI 58 [ prephitmp.771 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 75 74 76 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:148 (set (mem/s:DI (plus:DI (reg/f:DI 74 [ this ])
                (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64])
        (reg:DI 83)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (reg/v/f:DI 70 [ __tmp ])
            (reg:DI 58 [ prephitmp.771 ]))
        (nil)))

(jump_insn 76 75 77 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:148 (set (pc)
        (label_ref 122)) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

(barrier 77 76 78)

;; Start of basic block 7, registers live: (nil)
(code_label 78 77 79 7 39 "" [1 uses])

(note 79 78 80 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 80 79 81 7 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 402)

(insn 81 80 82 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (set (reg:DI 85)
        (mem/s/f:DI (plus:DI (reg/f:DI 74 [ this ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 82 81 83 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (parallel [
            (set (reg:DI 84)
                (minus:DI (reg:DI 85)
                    (reg:DI 64 [ D.38363 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (minus:DI (mem/s/f:DI (plus:DI (reg/f:DI 74 [ this ])
                    (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
            (reg:DI 64 [ D.38363 ]))
        (nil)))

(insn 83 82 84 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (parallel [
            (set (reg:DI 86)
                (ashiftrt:DI (reg:DI 84)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 84)
            (const_int 8 [0x8]))
        (nil)))

(insn 84 83 85 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:402 (set (reg:DI 72 [ D.35213 ])
        (reg:DI 86)) -1 (nil)
    (nil))

(note 85 84 86 7 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 150)

(insn 86 85 87 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:150 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 71 [ __xlen ])
            (reg:DI 72 [ D.35213 ]))) -1 (nil)
    (nil))

(jump_insn 87 86 89 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:150 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 98)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(note 89 87 90 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 90 89 91 8 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h") 300)

(insn 91 90 92 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 1 dx)
        (reg:DI 66 [ D.38333 ])) -1 (nil)
    (nil))

(insn 92 91 93 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 4 si)
        (reg:DI 65 [ D.38311 ])) -1 (nil)
    (nil))

(insn 93 92 94 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 5 di)
        (reg:DI 62 [ D.38345 ])) -1 (nil)
    (nil))

(call_insn 94 93 95 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memmove") [flags 0x41] <function_decl 0x2acebfb0b800 memmove>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 95 94 96 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 58 [ prephitmp.771 ])
                (ashift:DI (reg/v:DI 71 [ __xlen ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 96 95 97 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (pc)
        (label_ref 122)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

(barrier 97 96 98)

;; Start of basic block 9, registers live: (nil)
(code_label 98 97 99 9 46 "" [1 uses])

(note 99 98 100 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 100 99 101 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 87)
                (ashift:DI (reg:DI 72 [ D.35213 ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 101 100 102 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 88)
                (plus:DI (reg:DI 65 [ D.38311 ])
                    (reg:DI 87)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 102 101 103 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 89)
                (minus:DI (reg:DI 88)
                    (reg:DI 67 [ D.38332 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 103 102 104 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 1 dx)
        (reg:DI 89)) -1 (nil)
    (nil))

(insn 104 103 105 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 4 si)
        (reg:DI 65 [ D.38311 ])) -1 (nil)
    (nil))

(insn 105 104 106 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 5 di)
        (reg:DI 62 [ D.38345 ])) -1 (nil)
    (nil))

(call_insn 106 105 107 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memmove") [flags 0x41] <function_decl 0x2acebfb0b800 memmove>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(note 107 106 108 9 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 159)

(insn 108 107 109 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:159 (set (reg/v/f:DI 59 [ __result ])
        (mem/s/f:DI (plus:DI (reg/f:DI 74 [ this ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(note 109 108 110 9 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 704)

(insn 110 109 111 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:704 (parallel [
            (set (reg:DI 90)
                (minus:DI (reg/v/f:DI 59 [ __result ])
                    (mem/s/f:DI (reg/f:DI 74 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 111 110 112 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:704 (parallel [
            (set (reg:DI 92)
                (ashiftrt:DI (reg:DI 90)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 90)
            (const_int 8 [0x8]))
        (nil)))

(insn 112 111 113 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:704 (parallel [
            (set (reg:DI 93)
                (ashift:DI (reg:DI 92)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 113 112 114 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:704 (parallel [
            (set (reg:DI 68 [ D.36742 ])
                (plus:DI (reg:DI 93)
                    (mem/s/f:DI (reg/v/f:DI 75 [ __x ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 114 113 115 9 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h") 300)

(insn 115 114 116 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 95)
        (mem/s/f:DI (plus:DI (reg/v/f:DI 75 [ __x ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 116 115 117 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 94)
                (minus:DI (reg:DI 95)
                    (reg:DI 68 [ D.36742 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (minus:DI (mem/s/f:DI (plus:DI (reg/v/f:DI 75 [ __x ])
                    (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
            (reg:DI 68 [ D.36742 ]))
        (nil)))

(insn 117 116 118 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 1 dx)
        (reg:DI 94)) -1 (nil)
    (nil))

(insn 118 117 119 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 4 si)
        (reg:DI 68 [ D.36742 ])) -1 (nil)
    (nil))

(insn 119 118 120 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 5 di)
        (reg/v/f:DI 59 [ __result ])) -1 (nil)
    (nil))

(call_insn 120 119 121 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memmove") [flags 0x41] <function_decl 0x2acebfb0b800 memmove>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 121 120 122 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_algobase.h:300 (parallel [
            (set (reg:DI 58 [ prephitmp.771 ])
                (ashift:DI (reg/v:DI 71 [ __xlen ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(code_label 122 121 123 10 45 "" [2 uses])

(note 123 122 124 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 124 123 125 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 163)

(insn 125 124 126 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:163 (parallel [
            (set (reg:DI 96)
                (plus:DI (mem/s/f:DI (reg/f:DI 74 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])
                    (reg:DI 58 [ prephitmp.771 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 126 125 127 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:163 (set (mem/s:DI (plus:DI (reg/f:DI 74 [ this ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
        (reg:DI 96)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (mem/s/f:DI (reg/f:DI 74 [ this ]) [4 <variable>.D.34734._M_impl._M_start+0 S8 A64])
            (reg:DI 58 [ prephitmp.771 ]))
        (nil)))
;; End of basic block 10, registers live:
 (nil)

;; Start of basic block 11, registers live: (nil)
(code_label 127 126 128 11 37 "" [1 uses])

(note 128 127 129 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 129 128 132 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:163 (set (reg:DI 73 [ <result> ])
        (reg/f:DI 74 [ this ])) -1 (nil)
    (nil))

(note 132 129 133 11 NOTE_INSN_FUNCTION_END)

(note 133 132 135 11 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc") 165)

(insn 135 133 141 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:165 (set (reg/i:DI 0 ax)
        (reg:DI 73 [ <result> ])) -1 (nil)
    (nil))

(insn 141 135 0 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/vector.tcc:165 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)


;; Function void getAllKmers(int, std::vector<char*, std::allocator<char*> >&) (_Z11getAllKmersiRSt6vectorIPcSaIS0_EE)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Deleted label in block 7.
Deleted label in block 8.
Deleted label in block 12.
Deleted label in block 13.
Deleted label in block 17.
Deleted label in block 18.
Deleted label in block 22.
Forwarding edge 22->23 to 47 failed.
Forwarding edge 25->26 to 36 failed.
Deleted label in block 31.
Deleted label in block 32.
Forwarding edge 36->37 to 42 failed.
Deleted label in block 39.
Forwarding edge 42->43 to 76 failed.
Redirecting jump 329 from 50 to 76.
Deleted label in block 45.
Edge 47->50 redirected to 76
Deleted label in block 48.
Redirecting jump 368 from 50 to 76.
Merged 49 and 50 without moving.
Redirecting jump 376 from 75 to 76.
Deleted label in block 53.
Deleted label in block 54.
Forwarding edge 54->55 to 59 failed.
Deleted label in block 60.
Deleted label in block 61.
Forwarding edge 61->62 to 66 failed.
Deleted label in block 67.
Deleted label in block 68.
Forwarding edge 68->69 to 73 failed.
Forwarding edge 73->74 to 30 failed.
Deleting block 75.


try_optimize_cfg iteration 2

Forwarding edge 22->23 to 47 failed.
Forwarding edge 25->26 to 36 failed.
Forwarding edge 36->37 to 42 failed.
Forwarding edge 42->43 to 76 failed.
Forwarding edge 54->55 to 59 failed.
Forwarding edge 61->62 to 66 failed.
Forwarding edge 68->69 to 73 failed.
Forwarding edge 73->74 to 30 failed.


try_optimize_cfg iteration 1

Forwarding edge 21->22 to 48 failed.
Forwarding edge 24->25 to 35 failed.
Forwarding edge 35->36 to 42 failed.
Merged 37 and 38 without moving.
Forwarding edge 42->43 to 75 failed.
Forwarding edge 54->55 to 59 failed.
Forwarding edge 61->62 to 66 failed.
Forwarding edge 68->69 to 73 failed.
Forwarding edge 73->74 to 29 failed.


try_optimize_cfg iteration 2

Forwarding edge 21->22 to 48 failed.
Forwarding edge 24->25 to 35 failed.
Forwarding edge 35->36 to 42 failed.
Forwarding edge 42->43 to 75 failed.
Forwarding edge 54->55 to 59 failed.
Forwarding edge 61->62 to 66 failed.
Forwarding edge 68->69 to 73 failed.
Forwarding edge 73->74 to 29 failed.
(note 1 0 10 ("./DNAsequence.cc") 31)

;; Start of basic block 0, registers live: (nil)
(note 10 1 6 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 6 10 7 0 ./DNAsequence.cc:31 (set (reg/v:SI 108 [ k ])
        (reg:SI 5 di [ k ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./DNAsequence.cc:31 (set (reg/v/f:DI 109 [ kmers ])
        (reg:DI 4 si [ kmers ])) -1 (nil)
    (nil))

(note 8 7 12 0 NOTE_INSN_FUNCTION_BEG)

(note 12 8 13 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 87)

(insn 13 12 14 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:87 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 14 13 15 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:87 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 15 14 16 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:87 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -48 [0xffffffffffffffd0])) [4 nucs.D.34734._M_impl._M_end_of_storage+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 16 15 17 0 ("./DNAsequence.cc") 34)

(insn 17 16 18 0 ./DNAsequence.cc:34 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC1") [flags 0x2] <string_cst 0x2acec290f120>)) -1 (nil)
    (nil))

(call_insn 18 17 19 0 ./DNAsequence.cc:34 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 19 18 20 0 ./DNAsequence.cc:34 (set (reg/f:DI 110)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 110)
        (nil)))

(insn 20 19 21 0 ./DNAsequence.cc:34 (set (reg:DI 100 [ D.39277 ])
        (reg/f:DI 110)) -1 (nil)
    (nil))

(insn 21 20 22 0 ./DNAsequence.cc:34 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -32 [0xffffffffffffffe0])) [5 D.34747+0 S8 A64])
        (reg:DI 100 [ D.39277 ])) -1 (nil)
    (nil))

(note 22 21 23 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 23 22 24 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 101 [ D.39282 ])
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 24 23 25 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 101 [ D.39282 ])
            (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])) [4 nucs.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 25 24 27 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 42)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 913 [0x391])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 27 25 28 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 28 27 29 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg/v/f:DI 99 [ __p ])
        (reg:DI 101 [ D.39282 ])) -1 (nil)
    (nil))

(note 29 28 30 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 30 29 31 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 99 [ __p ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 31 30 33 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 35)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 33 31 34 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 34 33 35 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 99 [ __p ]) [5 S8 A64])
        (reg:DI 100 [ D.39277 ])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 35 34 36 3 54 "" [1 uses])

(note 36 35 37 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 37 36 38 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 38 37 39 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (reg:DI 111)
                (plus:DI (reg:DI 101 [ D.39282 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 39 38 40 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])
        (reg:DI 111)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (reg:DI 101 [ D.39282 ])
            (const_int 8 [0x8]))
        (nil)))

(jump_insn 40 39 41 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 52)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 41 40 42)

;; Start of basic block 4, registers live: (nil)
(code_label 42 41 43 4 52 "" [1 uses])

(note 43 42 44 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 44 43 45 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 45 44 46 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 98 [ D.39286 ])
        (reg:DI 101 [ D.39282 ])) -1 (nil)
    (nil))

(insn 46 45 47 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 112)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -32 [0xffffffffffffffe0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 47 46 48 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 113)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 48 47 49 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 112)) -1 (nil)
    (nil))

(insn 49 48 50 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 98 [ D.39286 ])) -1 (nil)
    (nil))

(insn 50 49 51 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg:DI 113)) -1 (nil)
    (nil))

(call_insn 51 50 52 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 3 [0x3])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 52 51 53 5 56 "" [1 uses])

(note 53 52 54 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 54 53 55 5 ("./DNAsequence.cc") 35)

(insn 55 54 56 5 ./DNAsequence.cc:35 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC2") [flags 0x2] <string_cst 0x2acec290f210>)) -1 (nil)
    (nil))

(call_insn 56 55 57 5 ./DNAsequence.cc:35 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 57 56 58 5 ./DNAsequence.cc:35 (set (reg/f:DI 114)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 114)
        (nil)))

(insn 58 57 59 5 ./DNAsequence.cc:35 (set (reg:DI 96 [ D.39307 ])
        (reg/f:DI 114)) -1 (nil)
    (nil))

(insn 59 58 60 5 ./DNAsequence.cc:35 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -24 [0xffffffffffffffe8])) [5 D.34748+0 S8 A64])
        (reg:DI 96 [ D.39307 ])) -1 (nil)
    (nil))

(note 60 59 61 5 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 61 60 62 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 97 [ D.39312 ])
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 62 61 63 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 97 [ D.39312 ])
            (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])) [4 nucs.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 63 62 65 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 80)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 913 [0x391])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(note 65 63 66 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 66 65 67 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg/v/f:DI 95 [ __p ])
        (reg:DI 97 [ D.39312 ])) -1 (nil)
    (nil))

(note 67 66 68 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 68 67 69 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 95 [ __p ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 69 68 71 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 73)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 71 69 72 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 72 71 73 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 95 [ __p ]) [5 S8 A64])
        (reg:DI 96 [ D.39307 ])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 73 72 74 8 59 "" [1 uses])

(note 74 73 75 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 75 74 76 8 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 76 75 77 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (reg:DI 115)
                (plus:DI (reg:DI 97 [ D.39312 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 77 76 78 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])
        (reg:DI 115)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (reg:DI 97 [ D.39312 ])
            (const_int 8 [0x8]))
        (nil)))

(jump_insn 78 77 79 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 90)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

(barrier 79 78 80)

;; Start of basic block 9, registers live: (nil)
(code_label 80 79 81 9 57 "" [1 uses])

(note 81 80 82 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(note 82 81 83 9 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 83 82 84 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 94 [ D.39316 ])
        (reg:DI 97 [ D.39312 ])) -1 (nil)
    (nil))

(insn 84 83 85 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 116)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 85 84 86 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 117)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 86 85 87 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 116)) -1 (nil)
    (nil))

(insn 87 86 88 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 94 [ D.39316 ])) -1 (nil)
    (nil))

(insn 88 87 89 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg:DI 117)) -1 (nil)
    (nil))

(call_insn 89 88 90 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 3 [0x3])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(code_label 90 89 91 10 61 "" [1 uses])

(note 91 90 92 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 92 91 93 10 ("./DNAsequence.cc") 36)

(insn 93 92 94 10 ./DNAsequence.cc:36 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2acec290f240>)) -1 (nil)
    (nil))

(call_insn 94 93 95 10 ./DNAsequence.cc:36 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 95 94 96 10 ./DNAsequence.cc:36 (set (reg/f:DI 118)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 118)
        (nil)))

(insn 96 95 97 10 ./DNAsequence.cc:36 (set (reg:DI 92 [ D.39337 ])
        (reg/f:DI 118)) -1 (nil)
    (nil))

(insn 97 96 98 10 ./DNAsequence.cc:36 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -16 [0xfffffffffffffff0])) [5 D.34749+0 S8 A64])
        (reg:DI 92 [ D.39337 ])) -1 (nil)
    (nil))

(note 98 97 99 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 99 98 100 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 93 [ D.39342 ])
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 100 99 101 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 93 [ D.39342 ])
            (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])) [4 nucs.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 101 100 103 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 118)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 913 [0x391])
        (nil)))
;; End of basic block 10, registers live:
 (nil)

;; Start of basic block 11, registers live: (nil)
(note 103 101 104 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 104 103 105 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg/v/f:DI 91 [ __p ])
        (reg:DI 93 [ D.39342 ])) -1 (nil)
    (nil))

(note 105 104 106 11 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 106 105 107 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 91 [ __p ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 107 106 109 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 111)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(note 109 107 110 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(insn 110 109 111 12 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 91 [ __p ]) [5 S8 A64])
        (reg:DI 92 [ D.39337 ])) -1 (nil)
    (nil))
;; End of basic block 12, registers live:
 (nil)

;; Start of basic block 13, registers live: (nil)
(code_label 111 110 112 13 64 "" [1 uses])

(note 112 111 113 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(note 113 112 114 13 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 114 113 115 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (reg:DI 119)
                (plus:DI (reg:DI 93 [ D.39342 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 115 114 116 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])
        (reg:DI 119)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (reg:DI 93 [ D.39342 ])
            (const_int 8 [0x8]))
        (nil)))

(jump_insn 116 115 117 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 128)) -1 (nil)
    (nil))
;; End of basic block 13, registers live:
 (nil)

(barrier 117 116 118)

;; Start of basic block 14, registers live: (nil)
(code_label 118 117 119 14 62 "" [1 uses])

(note 119 118 120 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(note 120 119 121 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 121 120 122 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 90 [ D.39346 ])
        (reg:DI 93 [ D.39342 ])) -1 (nil)
    (nil))

(insn 122 121 123 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 120)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -16 [0xfffffffffffffff0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 123 122 124 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 121)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 124 123 125 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 120)) -1 (nil)
    (nil))

(insn 125 124 126 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 90 [ D.39346 ])) -1 (nil)
    (nil))

(insn 126 125 127 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg:DI 121)) -1 (nil)
    (nil))

(call_insn 127 126 128 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 3 [0x3])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 14, registers live:
 (nil)

;; Start of basic block 15, registers live: (nil)
(code_label 128 127 129 15 66 "" [1 uses])

(note 129 128 130 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 130 129 131 15 ("./DNAsequence.cc") 37)

(insn 131 130 132 15 ./DNAsequence.cc:37 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2acec290f270>)) -1 (nil)
    (nil))

(call_insn 132 131 133 15 ./DNAsequence.cc:37 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 133 132 134 15 ./DNAsequence.cc:37 (set (reg/f:DI 122)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 122)
        (nil)))

(insn 134 133 135 15 ./DNAsequence.cc:37 (set (reg:DI 88 [ D.39367 ])
        (reg/f:DI 122)) -1 (nil)
    (nil))

(insn 135 134 136 15 ./DNAsequence.cc:37 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -8 [0xfffffffffffffff8])) [5 D.34750+0 S8 A64])
        (reg:DI 88 [ D.39367 ])) -1 (nil)
    (nil))

(note 136 135 137 15 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 137 136 138 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 89 [ D.39372 ])
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 138 137 139 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 89 [ D.39372 ])
            (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])) [4 nucs.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 139 138 141 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 156)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 913 [0x391])
        (nil)))
;; End of basic block 15, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(note 141 139 142 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(insn 142 141 143 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg/v/f:DI 87 [ __p ])
        (reg:DI 89 [ D.39372 ])) -1 (nil)
    (nil))

(note 143 142 144 16 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 144 143 145 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 87 [ __p ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 145 144 147 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 149)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(note 147 145 148 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(insn 148 147 149 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 87 [ __p ]) [5 S8 A64])
        (reg:DI 88 [ D.39367 ])) -1 (nil)
    (nil))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(code_label 149 148 150 18 69 "" [1 uses])

(note 150 149 151 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(note 151 150 152 18 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 152 151 153 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (reg:DI 123)
                (plus:DI (reg:DI 89 [ D.39372 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 153 152 154 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -56 [0xffffffffffffffc8])) [4 nucs.D.34734._M_impl._M_finish+0 S8 A64])
        (reg:DI 123)) -1 (nil)
    (expr_list:REG_EQUAL (plus:DI (reg:DI 89 [ D.39372 ])
            (const_int 8 [0x8]))
        (nil)))

(jump_insn 154 153 155 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 166)) -1 (nil)
    (nil))
;; End of basic block 18, registers live:
 (nil)

(barrier 155 154 156)

;; Start of basic block 19, registers live: (nil)
(code_label 156 155 157 19 67 "" [1 uses])

(note 157 156 158 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(note 158 157 159 19 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 159 158 160 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 86 [ D.39376 ])
        (reg:DI 89 [ D.39372 ])) -1 (nil)
    (nil))

(insn 160 159 161 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 124)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -8 [0xfffffffffffffff8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 161 160 162 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 125)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 162 161 163 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 124)) -1 (nil)
    (nil))

(insn 163 162 164 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 86 [ D.39376 ])) -1 (nil)
    (nil))

(insn 164 163 165 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg:DI 125)) -1 (nil)
    (nil))

(call_insn 165 164 166 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 3 [0x3])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 19, registers live:
 (nil)

;; Start of basic block 20, registers live: (nil)
(code_label 166 165 167 20 71 "" [1 uses])

(note 167 166 168 20 [bb 20] NOTE_INSN_BASIC_BLOCK)

(note 168 167 169 20 ("./DNAsequence.cc") 39)

(insn 169 168 170 20 ./DNAsequence.cc:39 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 108 [ k ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 170 169 172 20 ./DNAsequence.cc:39 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 180)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7100 [0x1bbc])
        (nil)))
;; End of basic block 20, registers live:
 (nil)

;; Start of basic block 21, registers live: (nil)
(note 172 170 173 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(note 173 172 174 21 ("./DNAsequence.cc") 41)

(insn 174 173 175 21 ./DNAsequence.cc:41 (parallel [
            (set (reg:DI 126)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 175 174 176 21 ./DNAsequence.cc:41 (set (reg:DI 4 si)
        (reg:DI 126)) -1 (nil)
    (nil))

(insn 176 175 177 21 ./DNAsequence.cc:41 (set (reg:DI 5 di)
        (reg/v/f:DI 109 [ kmers ])) -1 (nil)
    (nil))

(call_insn 177 176 605 21 ./DNAsequence.cc:41 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EEaSERKS2_") [flags 0x1] <function_decl 0x2acec22bcd00 operator=>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 3 [0x3])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 21, registers live:
 (nil)

;; Start of basic block 22, registers live: (nil)
(note 605 177 178 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(jump_insn 178 605 179 22 ./DNAsequence.cc:41 (set (pc)
        (label_ref 354)) -1 (nil)
    (nil))
;; End of basic block 22, registers live:
 (nil)

(barrier 179 178 180)

;; Start of basic block 23, registers live: (nil)
(code_label 180 179 181 23 72 "" [1 uses])

(note 181 180 182 23 [bb 23] NOTE_INSN_BASIC_BLOCK)

(note 182 181 183 23 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 87)

(insn 183 182 184 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:87 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 184 183 185 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:87 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -88 [0xffffffffffffffa8])) [4 km1mers.D.34734._M_impl._M_finish+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 185 184 186 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:87 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -80 [0xffffffffffffffb0])) [4 km1mers.D.34734._M_impl._M_end_of_storage+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 186 185 187 23 ("./DNAsequence.cc") 46)

(insn 187 186 188 23 ./DNAsequence.cc:46 (parallel [
            (set (reg:DI 127)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -96 [0xffffffffffffffa0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 188 187 189 23 ./DNAsequence.cc:46 (parallel [
            (set (reg:SI 128)
                (plus:SI (reg/v:SI 108 [ k ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 189 188 190 23 ./DNAsequence.cc:46 (set (reg:DI 4 si)
        (reg:DI 127)) -1 (nil)
    (nil))

(insn 190 189 191 23 ./DNAsequence.cc:46 (set (reg:SI 5 di)
        (reg:SI 128)) -1 (nil)
    (nil))

(call_insn 191 190 192 23 ./DNAsequence.cc:46 (call (mem:QI (symbol_ref:DI ("_Z11getAllKmersiRSt6vectorIPcSaIS0_EE") [flags 0x3] <function_decl 0x2acec1315200 getAllKmers>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 23, registers live:
 (nil)

;; Start of basic block 24, registers live: (nil)
(note 192 191 193 24 [bb 24] NOTE_INSN_BASIC_BLOCK)

(note 193 192 194 24 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 809)

(insn 194 193 195 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 130)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 195 194 196 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 131)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -88 [0xffffffffffffffa8])) [4 km1mers.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 196 195 197 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (parallel [
            (set (reg:DI 129)
                (minus:DI (reg:DI 131)
                    (reg:DI 130)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (minus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -88 [0xffffffffffffffa8])) [4 km1mers.D.34734._M_impl._M_finish+0 S8 A64])
            (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64]))
        (nil)))

(insn 197 196 198 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (parallel [
            (set (reg:DI 132)
                (ashiftrt:DI (reg:DI 129)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 129)
            (const_int 8 [0x8]))
        (nil)))

(insn 198 197 199 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 85 [ D.39411 ])
        (reg:DI 132)) -1 (nil)
    (nil))

(note 199 198 200 24 ("./DNAsequence.cc") 52)

(insn 200 199 201 24 ./DNAsequence.cc:52 (set (reg:CCNO 17 flags)
        (compare:CCNO (subreg:SI (reg:DI 85 [ D.39411 ]) 0)
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 201 200 204 24 ./DNAsequence.cc:52 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 378)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9667 [0x25c3])
        (nil)))
;; End of basic block 24, registers live:
 (nil)

;; Start of basic block 25, registers live: (nil)
(note 204 201 202 25 [bb 25] NOTE_INSN_BASIC_BLOCK)

(jump_insn 202 204 203 25 ./DNAsequence.cc:52 (set (pc)
        (label_ref 266)) -1 (nil)
    (nil))
;; End of basic block 25, registers live:
 (nil)

(barrier 203 202 205)

;; Start of basic block 26, registers live: (nil)
(code_label 205 203 206 26 77 "" [4 uses])

(note 206 205 207 26 [bb 26] NOTE_INSN_BASIC_BLOCK)

(note 207 206 208 26 ("./DNAsequence.cc") 59)

(insn 208 207 209 26 ./DNAsequence.cc:59 (set (reg:DI 2 cx)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2acec0dfe370 stderr>) [20 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 209 208 210 26 ./DNAsequence.cc:59 (set (reg:DI 1 dx)
        (const_int 14 [0xe])) -1 (nil)
    (nil))

(insn 210 209 211 26 ./DNAsequence.cc:59 (set (reg:DI 4 si)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(insn 211 210 212 26 ./DNAsequence.cc:59 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC5") [flags 0x2] <string_cst 0x2acec2999740>)) -1 (nil)
    (nil))

(call_insn 212 211 213 26 ./DNAsequence.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fwrite") [flags 0x41] <function_decl 0x2acebfb1d000 __builtin_fwrite>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 26, registers live:
 (nil)

;; Start of basic block 27, registers live: (nil)
(note 213 212 214 27 [bb 27] NOTE_INSN_BASIC_BLOCK)

(insn 214 213 215 27 ./DNAsequence.cc:59 (set (reg:SI 2 cx)
        (const_int 59 [0x3b])) -1 (nil)
    (nil))

(insn 215 214 216 27 ./DNAsequence.cc:59 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC6") [flags 0x2] <string_cst 0x2acec2454880>)) -1 (nil)
    (nil))

(insn 216 215 217 27 ./DNAsequence.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC7") [flags 0x2] <string_cst 0x2acec2999880>)) -1 (nil)
    (nil))

(insn 217 216 218 27 ./DNAsequence.cc:59 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2acec0dfe370 stderr>) [20 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 218 217 219 27 ./DNAsequence.cc:59 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 219 218 220 27 ./DNAsequence.cc:59 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2acebfb17d00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))
;; End of basic block 27, registers live:
 (nil)

;; Start of basic block 28, registers live: (nil)
(note 220 219 221 28 [bb 28] NOTE_INSN_BASIC_BLOCK)

(insn 221 220 222 28 ./DNAsequence.cc:59 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 222 221 223 28 ./DNAsequence.cc:59 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2acebfb3bb00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 28, registers live:
 (nil)

(barrier 223 222 224)

;; Start of basic block 29, registers live: (nil)
(code_label 224 223 225 29 78 "" [1 uses])

(note 225 224 226 29 [bb 29] NOTE_INSN_BASIC_BLOCK)

(note 226 225 227 29 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 227 226 228 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 84 [ D.39438 ])
        (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 228 227 229 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 84 [ D.39438 ])
            (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                    (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 229 228 231 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 245)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 913 [0x391])
        (nil)))
;; End of basic block 29, registers live:
 (nil)

;; Start of basic block 30, registers live: (nil)
(note 231 229 232 30 [bb 30] NOTE_INSN_BASIC_BLOCK)

(insn 232 231 233 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg/v/f:DI 83 [ __p ])
        (reg:DI 84 [ D.39438 ])) -1 (nil)
    (nil))

(note 233 232 234 30 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 234 233 235 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 83 [ __p ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 235 234 237 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 239)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1900 [0x76c])
        (nil)))
;; End of basic block 30, registers live:
 (nil)

;; Start of basic block 31, registers live: (nil)
(note 237 235 238 31 [bb 31] NOTE_INSN_BASIC_BLOCK)

(insn 238 237 239 31 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 83 [ __p ]) [5 S8 A64])
        (reg:DI 63 [ kmer.1017 ])) -1 (nil)
    (nil))
;; End of basic block 31, registers live:
 (nil)

;; Start of basic block 32, registers live: (nil)
(code_label 239 238 240 32 81 "" [1 uses])

(note 240 239 241 32 [bb 32] NOTE_INSN_BASIC_BLOCK)

(note 241 240 242 32 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 242 241 243 32 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (mem/s:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                        (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                (plus:DI (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                            (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 243 242 244 32 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 254)) -1 (nil)
    (nil))
;; End of basic block 32, registers live:
 (nil)

(barrier 244 243 245)

;; Start of basic block 33, registers live: (nil)
(code_label 245 244 246 33 79 "" [1 uses])

(note 246 245 247 33 [bb 33] NOTE_INSN_BASIC_BLOCK)

(note 247 246 248 33 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 248 247 249 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 82 [ D.39442 ])
        (reg:DI 84 [ D.39438 ])) -1 (nil)
    (nil))

(insn 249 248 250 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 133)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -40 [0xffffffffffffffd8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 250 249 251 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 133)) -1 (nil)
    (nil))

(insn 251 250 252 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 82 [ D.39442 ])) -1 (nil)
    (nil))

(insn 252 251 253 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg/v/f:DI 109 [ kmers ])) -1 (nil)
    (nil))

(call_insn 253 252 254 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 33, registers live:
 (nil)

;; Start of basic block 34, registers live: (nil)
(code_label 254 253 255 34 83 "" [1 uses])

(note 255 254 256 34 [bb 34] NOTE_INSN_BASIC_BLOCK)

(note 256 255 257 34 ("./DNAsequence.cc") 62)

(insn 257 256 258 34 ./DNAsequence.cc:62 (set (reg:DI 134)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 258 257 259 34 ./DNAsequence.cc:62 (set (reg:DI 135)
        (mem/f:DI (plus:DI (reg:DI 71 [ pretmp.977 ])
                (reg:DI 134)) [5 S8 A64])) -1 (nil)
    (nil))

(insn 259 258 260 34 ./DNAsequence.cc:62 (set (reg:DI 5 di)
        (reg:DI 135)) -1 (nil)
    (nil))

(call_insn 260 259 261 34 ./DNAsequence.cc:62 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2acebff7a900 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(note 261 260 262 34 ("./DNAsequence.cc") 52)

(insn 262 261 263 34 ./DNAsequence.cc:52 (parallel [
            (set (reg/v:SI 106 [ i ])
                (plus:SI (reg/v:SI 106 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 263 262 264 34 ./DNAsequence.cc:52 (parallel [
            (set (reg:DI 70 [ ivtmp.990 ])
                (plus:DI (reg:DI 70 [ ivtmp.990 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 264 263 265 34 ./DNAsequence.cc:52 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 106 [ i ])
            (subreg:SI (reg:DI 85 [ D.39411 ]) 0))) -1 (nil)
    (nil))

(jump_insn 265 264 266 34 ./DNAsequence.cc:52 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 382)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9667 [0x25c3])
        (nil)))
;; End of basic block 34, registers live:
 (nil)

;; Start of basic block 35, registers live: (nil)
(code_label 266 265 267 35 76 "" [1 uses])

(note 267 266 268 35 [bb 35] NOTE_INSN_BASIC_BLOCK)

(note 268 267 269 35 ("./DNAsequence.cc") 66)

(insn 269 268 270 35 ./DNAsequence.cc:66 (set (reg/f:DI 136)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 270 269 271 35 ./DNAsequence.cc:66 (set (reg:DI 5 di)
        (mem/f:DI (reg/f:DI 136) [5 S8 A64])) -1 (nil)
    (nil))

(call_insn 271 270 272 35 ./DNAsequence.cc:66 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2acebff7a900 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 272 271 273 35 ./DNAsequence.cc:66 (set (reg/f:DI 137)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 273 272 274 35 ./DNAsequence.cc:66 (set (reg:DI 138)
        (mem/f:DI (plus:DI (reg/f:DI 137)
                (const_int 8 [0x8])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 274 273 275 35 ./DNAsequence.cc:66 (set (reg:DI 5 di)
        (reg:DI 138)) -1 (nil)
    (nil))

(call_insn 275 274 276 35 ./DNAsequence.cc:66 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2acebff7a900 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 276 275 277 35 ./DNAsequence.cc:66 (set (reg/f:DI 139)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 277 276 278 35 ./DNAsequence.cc:66 (set (reg:DI 140)
        (mem/f:DI (plus:DI (reg/f:DI 139)
                (const_int 16 [0x10])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 278 277 279 35 ./DNAsequence.cc:66 (set (reg:DI 5 di)
        (reg:DI 140)) -1 (nil)
    (nil))

(call_insn 279 278 280 35 ./DNAsequence.cc:66 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2acebff7a900 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 280 279 281 35 ./DNAsequence.cc:66 (set (reg/f:DI 141)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 281 280 282 35 ./DNAsequence.cc:66 (set (reg:DI 142)
        (mem/f:DI (plus:DI (reg/f:DI 141)
                (const_int 24 [0x18])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 282 281 283 35 ./DNAsequence.cc:66 (set (reg:DI 5 di)
        (reg:DI 142)) -1 (nil)
    (nil))

(call_insn 283 282 284 35 ./DNAsequence.cc:66 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2acebff7a900 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(note 284 283 285 35 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 272)

(insn 285 284 286 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:272 (parallel [
            (set (reg/f:DI 79 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -96 [0xffffffffffffffa0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 286 285 287 35 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 119)

(insn 287 286 288 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:119 (set (reg:DI 78 [ D.39569 ])
        (mem/s/f:DI (reg/f:DI 79 [ this ]) [4 <variable>._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 288 287 289 35 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 132)

(insn 289 288 290 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 78 [ D.39569 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 290 289 293 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 316)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 6463 [0x193f])
        (nil)))
;; End of basic block 35, registers live:
 (nil)

;; Start of basic block 36, registers live: (nil)
(note 293 290 291 36 [bb 36] NOTE_INSN_BASIC_BLOCK)

(jump_insn 291 293 292 36 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (label_ref 320)) -1 (nil)
    (nil))
;; End of basic block 36, registers live:
 (nil)

(barrier 292 291 616)

;; Start of basic block 37, registers live: (nil)
(code_label/s 616 292 619 37 115 "" [1 uses])

(note 619 616 617 37 [bb 37] NOTE_INSN_BASIC_BLOCK)

(insn 617 619 618 37 (set (reg:DI 144)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 618 617 294 37 (set (reg:DI 143)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 294 618 296 37 "" NOTE_INSN_DELETED_LABEL 87)

(insn 296 294 297 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:SI 104 [ save_filt.122 ])
        (subreg:SI (reg:DI 143) 0)) -1 (nil)
    (nil))

(insn 297 296 298 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:DI 105 [ save_eptr.121 ])
        (reg:DI 144)) -1 (nil)
    (nil))

(note 298 297 299 37 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 272)

(insn 299 298 300 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:272 (parallel [
            (set (reg/f:DI 81 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -96 [0xffffffffffffffa0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 300 299 301 37 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 119)

(insn 301 300 302 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:119 (set (reg:DI 80 [ D.39498 ])
        (mem/s/f:DI (reg/f:DI 81 [ this ]) [4 <variable>._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 302 301 303 37 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 132)

(insn 303 302 304 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 80 [ D.39498 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 304 303 306 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 310)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3537 [0xdd1])
        (nil)))
;; End of basic block 37, registers live:
 (nil)

;; Start of basic block 39, registers live: (nil)
(note 306 304 307 39 [bb 39] NOTE_INSN_BASIC_BLOCK)

(note 307 306 308 39 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 94)

(insn 308 307 309 39 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg:DI 80 [ D.39498 ])) -1 (nil)
    (nil))

(call_insn 309 308 310 39 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_ZdlPv") [flags 0x41] <function_decl 0x2acebfb9a900 operator delete>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 39, registers live:
 (nil)

;; Start of basic block 40, registers live: (nil)
(code_label 310 309 311 40 88 "" [1 uses])

(note 311 310 312 40 [bb 40] NOTE_INSN_BASIC_BLOCK)

(insn 312 311 313 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 144)
        (reg:DI 105 [ save_eptr.121 ])) -1 (nil)
    (nil))

(insn 313 312 611 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 143)
        (sign_extend:DI (reg:SI 104 [ save_filt.122 ]))) -1 (nil)
    (nil))

(jump_insn 611 313 315 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (pc)
        (label_ref 332)) -1 (nil)
    (nil))
;; End of basic block 40, registers live:
 (nil)

(barrier 315 611 316)

;; Start of basic block 41, registers live: (nil)
(code_label 316 315 317 41 85 "" [1 uses])

(note 317 316 318 41 [bb 41] NOTE_INSN_BASIC_BLOCK)

(insn 318 317 319 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg:DI 78 [ D.39569 ])) -1 (nil)
    (nil))

(call_insn 319 318 320 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_ZdlPv") [flags 0x41] <function_decl 0x2acebfb9a900 operator delete>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 41, registers live:
 (nil)

;; Start of basic block 42, registers live: (nil)
(code_label 320 319 321 42 86 "" [1 uses])

(note 321 320 322 42 [bb 42] NOTE_INSN_BASIC_BLOCK)

(note 322 321 323 42 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 272)

(insn 323 322 324 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:272 (parallel [
            (set (reg/f:DI 73 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 324 323 325 42 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 119)

(insn 325 324 326 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:119 (set (reg:DI 72 [ D.39782 ])
        (mem/s/f:DI (reg/f:DI 73 [ this ]) [4 <variable>._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 326 325 327 42 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 132)

(insn 327 326 328 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 72 [ D.39782 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 328 327 331 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 370)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7876 [0x1ec4])
        (nil)))
;; End of basic block 42, registers live:
 (nil)

;; Start of basic block 43, registers live: (nil)
(note 331 328 329 43 [bb 43] NOTE_INSN_BASIC_BLOCK)

(jump_insn 329 331 330 43 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (label_ref:DI 603)) 548 {jump} (nil)
    (nil))
;; End of basic block 43, registers live:
 (nil)

(barrier 330 329 620)

;; Start of basic block 44, registers live: (nil)
(code_label/s 620 330 623 44 116 "" [1 uses])

(note 623 620 621 44 [bb 44] NOTE_INSN_BASIC_BLOCK)

(insn 621 623 622 44 (set (reg:DI 144)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 622 621 332 44 (set (reg:DI 143)
        (reg:DI 1 dx)) -1 (nil)
    (nil))
;; End of basic block 44, registers live:
 (nil)

;; Start of basic block 45, registers live: (nil)
(code_label/s 332 622 333 45 92 "" [2 uses])

(note 333 332 334 45 [bb 45] NOTE_INSN_BASIC_BLOCK)

(insn 334 333 335 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:SI 102 [ save_filt.124 ])
        (subreg:SI (reg:DI 143) 0)) -1 (nil)
    (nil))

(insn 335 334 336 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:DI 103 [ save_eptr.123 ])
        (reg:DI 144)) -1 (nil)
    (nil))

(note 336 335 337 45 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 272)

(insn 337 336 338 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:272 (parallel [
            (set (reg/f:DI 77 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 338 337 339 45 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 119)

(insn 339 338 340 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:119 (set (reg:DI 76 [ D.39640 ])
        (mem/s/f:DI (reg/f:DI 77 [ this ]) [4 <variable>._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 340 339 341 45 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 132)

(insn 341 340 342 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 76 [ D.39640 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 342 341 344 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 348)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3537 [0xdd1])
        (nil)))
;; End of basic block 45, registers live:
 (nil)

;; Start of basic block 46, registers live: (nil)
(note 344 342 345 46 [bb 46] NOTE_INSN_BASIC_BLOCK)

(note 345 344 346 46 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 94)

(insn 346 345 347 46 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg:DI 76 [ D.39640 ])) -1 (nil)
    (nil))

(call_insn 347 346 348 46 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_ZdlPv") [flags 0x41] <function_decl 0x2acebfb9a900 operator delete>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 46, registers live:
 (nil)

;; Start of basic block 47, registers live: (nil)
(code_label 348 347 349 47 93 "" [1 uses])

(note 349 348 350 47 [bb 47] NOTE_INSN_BASIC_BLOCK)

(insn 350 349 351 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 144)
        (reg:DI 103 [ save_eptr.123 ])) -1 (nil)
    (nil))

(insn 351 350 613 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 143)
        (sign_extend:DI (reg:SI 102 [ save_filt.124 ]))) -1 (nil)
    (nil))

(insn 613 351 614 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg:DI 144)) -1 (nil)
    (nil))

(call_insn 614 613 353 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 47, registers live:
 (nil)

(barrier 353 614 354)

;; Start of basic block 48, registers live: (nil)
(code_label 354 353 355 48 74 "" [1 uses])

(note 355 354 356 48 [bb 48] NOTE_INSN_BASIC_BLOCK)

(note 356 355 357 48 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 272)

(insn 357 356 358 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:272 (parallel [
            (set (reg/f:DI 75 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 358 357 359 48 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 119)

(insn 359 358 360 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:119 (set (reg:DI 74 [ D.39711 ])
        (mem/s/f:DI (reg/f:DI 75 [ this ]) [4 <variable>._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(note 360 359 361 48 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 132)

(insn 361 360 362 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 74 [ D.39711 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 362 361 364 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:132 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 603)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 2124 [0x84c])
        (nil)))
;; End of basic block 48, registers live:
 (nil)

;; Start of basic block 49, registers live: (nil)
(note 364 362 365 49 [bb 49] NOTE_INSN_BASIC_BLOCK)

(note 365 364 366 49 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 94)

(insn 366 365 367 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg:DI 74 [ D.39711 ])) -1 (nil)
    (nil))

(call_insn 367 366 368 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_ZdlPv") [flags 0x41] <function_decl 0x2acebfb9a900 operator delete>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 368 367 369 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (pc)
        (label_ref:DI 603)) 548 {jump} (nil)
    (nil))
;; End of basic block 49, registers live:
 (nil)

(barrier 369 368 370)

;; Start of basic block 50, registers live: (nil)
(code_label 370 369 371 50 90 "" [1 uses])

(note 371 370 372 50 [bb 50] NOTE_INSN_BASIC_BLOCK)

(insn 372 371 373 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 5 di)
        (reg:DI 72 [ D.39782 ])) -1 (nil)
    (nil))

(call_insn 373 372 376 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (call (mem:QI (symbol_ref:DI ("_ZdlPv") [flags 0x41] <function_decl 0x2acebfb9a900 operator delete>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 376 373 377 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (pc)
        (label_ref:DI 603)) 548 {jump} (nil)
    (nil))
;; End of basic block 50, registers live:
 (nil)

(barrier 377 376 378)

;; Start of basic block 51, registers live: (nil)
(code_label 378 377 379 51 75 "" [1 uses])

(note 379 378 380 51 [bb 51] NOTE_INSN_BASIC_BLOCK)

(insn 380 379 381 51 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg/v:SI 106 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 381 380 382 51 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 70 [ ivtmp.990 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 51, registers live:
 (nil)

;; Start of basic block 52, registers live: (nil)
(code_label 382 381 383 52 84 "" [1 uses])

(note 383 382 384 52 [bb 52] NOTE_INSN_BASIC_BLOCK)

(insn 384 383 385 52 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:94 (set (reg:DI 71 [ pretmp.977 ])
        (reg:DI 70 [ ivtmp.990 ])) -1 (nil)
    (nil))

(note 385 384 386 52 ("./DNAsequence.cc") 56)

(insn 386 385 387 52 ./DNAsequence.cc:56 (set (mem/s:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -208 [0xffffffffffffff30])) [0 str+0 S1 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 387 386 388 52 ("./DNAsequence.cc") 57)

(insn 388 387 389 52 ./DNAsequence.cc:57 (set (reg/f:DI 145)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 389 388 390 52 ./DNAsequence.cc:57 (parallel [
            (set (reg:DI 146)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 390 389 391 52 ./DNAsequence.cc:57 (set (reg:DI 4 si)
        (mem/f:DI (reg/f:DI 145) [5 S8 A64])) -1 (nil)
    (nil))

(insn 391 390 392 52 ./DNAsequence.cc:57 (set (reg:DI 5 di)
        (reg:DI 146)) -1 (nil)
    (nil))

(call_insn 392 391 393 52 ./DNAsequence.cc:57 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 393 392 394 52 ("./DNAsequence.cc") 58)

(insn 394 393 395 52 ./DNAsequence.cc:58 (set (reg:DI 147)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 395 394 396 52 ./DNAsequence.cc:58 (set (reg:DI 148)
        (mem/f:DI (plus:DI (reg:DI 71 [ pretmp.977 ])
                (reg:DI 147)) [5 S8 A64])) -1 (nil)
    (nil))

(insn 396 395 397 52 ./DNAsequence.cc:58 (parallel [
            (set (reg:DI 149)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 397 396 398 52 ./DNAsequence.cc:58 (set (reg:DI 4 si)
        (reg:DI 148)) -1 (nil)
    (nil))

(insn 398 397 399 52 ./DNAsequence.cc:58 (set (reg:DI 5 di)
        (reg:DI 149)) -1 (nil)
    (nil))

(call_insn 399 398 400 52 ./DNAsequence.cc:58 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 400 399 401 52 ("./DNAsequence.cc") 59)

(insn 401 400 402 52 ./DNAsequence.cc:59 (parallel [
            (set (reg:DI 150)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 402 401 403 52 ./DNAsequence.cc:59 (set (reg:DI 5 di)
        (reg:DI 150)) -1 (nil)
    (nil))

(call_insn 403 402 404 52 ./DNAsequence.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 404 403 405 52 ./DNAsequence.cc:59 (set (reg/f:DI 151)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 151)
        (nil)))

(insn 405 404 406 52 ./DNAsequence.cc:59 (set (reg:DI 107 [ kmer.94 ])
        (reg/f:DI 151)) -1 (nil)
    (nil))

(insn 406 405 407 52 ./DNAsequence.cc:59 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -40 [0xffffffffffffffd8])) [5 kmer+0 S8 A64])
        (reg:DI 107 [ kmer.94 ])) -1 (nil)
    (nil))

(insn 407 406 408 52 ./DNAsequence.cc:59 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 107 [ kmer.94 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 408 407 410 52 ./DNAsequence.cc:59 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 205)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 52, registers live:
 (nil)

;; Start of basic block 53, registers live: (nil)
(note 410 408 411 53 [bb 53] NOTE_INSN_BASIC_BLOCK)

(note 411 410 412 53 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 412 411 413 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 66 [ temp.1014 ])
        (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 413 412 414 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 66 [ temp.1014 ])
            (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                    (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 414 413 416 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 426)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9087 [0x237f])
        (nil)))
;; End of basic block 53, registers live:
 (nil)

;; Start of basic block 54, registers live: (nil)
(note 416 414 417 54 [bb 54] NOTE_INSN_BASIC_BLOCK)

(note 417 416 418 54 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 418 417 419 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 82 [ D.39442 ])
        (reg:DI 66 [ temp.1014 ])) -1 (nil)
    (nil))

(insn 419 418 420 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 152)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -40 [0xffffffffffffffd8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 420 419 421 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 152)) -1 (nil)
    (nil))

(insn 421 420 422 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 82 [ D.39442 ])) -1 (nil)
    (nil))

(insn 422 421 423 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg/v/f:DI 109 [ kmers ])) -1 (nil)
    (nil))

(call_insn 423 422 606 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 54, registers live:
 (nil)

;; Start of basic block 55, registers live: (nil)
(note 606 423 424 55 [bb 55] NOTE_INSN_BASIC_BLOCK)

(jump_insn 424 606 425 55 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (pc)
        (label_ref 444)) -1 (nil)
    (nil))
;; End of basic block 55, registers live:
 (nil)

(barrier 425 424 426)

;; Start of basic block 56, registers live: (nil)
(code_label 426 425 427 56 97 "" [1 uses])

(note 427 426 428 56 [bb 56] NOTE_INSN_BASIC_BLOCK)

(insn 428 427 429 56 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg/v/f:DI 69 [ __p.1008 ])
        (reg:DI 66 [ temp.1014 ])) -1 (nil)
    (nil))

(note 429 428 430 56 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 430 429 431 56 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 69 [ __p.1008 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 431 430 432 56 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 438)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8100 [0x1fa4])
        (nil)))
;; End of basic block 56, registers live:
 (nil)

;; Start of basic block 57, registers live: (nil)
(code_label 432 431 433 57 101 "" [1 uses])

(note 433 432 434 57 [bb 57] NOTE_INSN_BASIC_BLOCK)

(note 434 433 435 57 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 435 434 436 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (mem/s:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                        (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                (plus:DI (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                            (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 436 435 437 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 444)) -1 (nil)
    (nil))
;; End of basic block 57, registers live:
 (nil)

(barrier 437 436 438)

;; Start of basic block 58, registers live: (nil)
(code_label 438 437 439 58 100 "" [1 uses])

(note 439 438 440 58 [bb 58] NOTE_INSN_BASIC_BLOCK)

(note 440 439 441 58 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 441 440 442 58 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 69 [ __p.1008 ]) [5 S8 A64])
        (reg:DI 107 [ kmer.94 ])) -1 (nil)
    (nil))

(jump_insn 442 441 443 58 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (label_ref 432)) -1 (nil)
    (nil))
;; End of basic block 58, registers live:
 (nil)

(barrier 443 442 444)

;; Start of basic block 59, registers live: (nil)
(code_label 444 443 445 59 99 "" [2 uses])

(note 445 444 446 59 [bb 59] NOTE_INSN_BASIC_BLOCK)

(note 446 445 447 59 ("./DNAsequence.cc") 56)

(insn 447 446 448 59 ./DNAsequence.cc:56 (set (mem/s:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -208 [0xffffffffffffff30])) [0 str+0 S1 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 448 447 449 59 ("./DNAsequence.cc") 57)

(insn 449 448 450 59 ./DNAsequence.cc:57 (set (reg/f:DI 153)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 450 449 451 59 ./DNAsequence.cc:57 (set (reg:DI 154)
        (mem/f:DI (plus:DI (reg/f:DI 153)
                (const_int 8 [0x8])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 451 450 452 59 ./DNAsequence.cc:57 (parallel [
            (set (reg:DI 155)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 452 451 453 59 ./DNAsequence.cc:57 (set (reg:DI 4 si)
        (reg:DI 154)) -1 (nil)
    (nil))

(insn 453 452 454 59 ./DNAsequence.cc:57 (set (reg:DI 5 di)
        (reg:DI 155)) -1 (nil)
    (nil))

(call_insn 454 453 455 59 ./DNAsequence.cc:57 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 455 454 456 59 ("./DNAsequence.cc") 58)

(insn 456 455 457 59 ./DNAsequence.cc:58 (set (reg:DI 156)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 457 456 458 59 ./DNAsequence.cc:58 (set (reg:DI 157)
        (mem/f:DI (plus:DI (reg:DI 71 [ pretmp.977 ])
                (reg:DI 156)) [5 S8 A64])) -1 (nil)
    (nil))

(insn 458 457 459 59 ./DNAsequence.cc:58 (parallel [
            (set (reg:DI 158)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 459 458 460 59 ./DNAsequence.cc:58 (set (reg:DI 4 si)
        (reg:DI 157)) -1 (nil)
    (nil))

(insn 460 459 461 59 ./DNAsequence.cc:58 (set (reg:DI 5 di)
        (reg:DI 158)) -1 (nil)
    (nil))

(call_insn 461 460 462 59 ./DNAsequence.cc:58 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 462 461 463 59 ("./DNAsequence.cc") 59)

(insn 463 462 464 59 ./DNAsequence.cc:59 (parallel [
            (set (reg:DI 159)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 464 463 465 59 ./DNAsequence.cc:59 (set (reg:DI 5 di)
        (reg:DI 159)) -1 (nil)
    (nil))

(call_insn 465 464 466 59 ./DNAsequence.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 466 465 467 59 ./DNAsequence.cc:59 (set (reg/f:DI 160)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 160)
        (nil)))

(insn 467 466 468 59 ./DNAsequence.cc:59 (set (reg:DI 62 [ kmer.1018 ])
        (reg/f:DI 160)) -1 (nil)
    (nil))

(insn 468 467 469 59 ./DNAsequence.cc:59 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -40 [0xffffffffffffffd8])) [5 kmer+0 S8 A64])
        (reg:DI 62 [ kmer.1018 ])) -1 (nil)
    (nil))

(insn 469 468 470 59 ./DNAsequence.cc:59 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 62 [ kmer.1018 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 470 469 472 59 ./DNAsequence.cc:59 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 205)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 59, registers live:
 (nil)

;; Start of basic block 60, registers live: (nil)
(note 472 470 473 60 [bb 60] NOTE_INSN_BASIC_BLOCK)

(note 473 472 474 60 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 474 473 475 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 65 [ temp.1015 ])
        (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 475 474 476 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 65 [ temp.1015 ])
            (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                    (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 476 475 478 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 488)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9087 [0x237f])
        (nil)))
;; End of basic block 60, registers live:
 (nil)

;; Start of basic block 61, registers live: (nil)
(note 478 476 479 61 [bb 61] NOTE_INSN_BASIC_BLOCK)

(note 479 478 480 61 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 480 479 481 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 82 [ D.39442 ])
        (reg:DI 65 [ temp.1015 ])) -1 (nil)
    (nil))

(insn 481 480 482 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 161)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -40 [0xffffffffffffffd8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 482 481 483 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 161)) -1 (nil)
    (nil))

(insn 483 482 484 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 82 [ D.39442 ])) -1 (nil)
    (nil))

(insn 484 483 485 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg/v/f:DI 109 [ kmers ])) -1 (nil)
    (nil))

(call_insn 485 484 607 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 61, registers live:
 (nil)

;; Start of basic block 62, registers live: (nil)
(note 607 485 486 62 [bb 62] NOTE_INSN_BASIC_BLOCK)

(jump_insn 486 607 487 62 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (pc)
        (label_ref 506)) -1 (nil)
    (nil))
;; End of basic block 62, registers live:
 (nil)

(barrier 487 486 488)

;; Start of basic block 63, registers live: (nil)
(code_label 488 487 489 63 103 "" [1 uses])

(note 489 488 490 63 [bb 63] NOTE_INSN_BASIC_BLOCK)

(insn 490 489 491 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg/v/f:DI 68 [ __p.1009 ])
        (reg:DI 65 [ temp.1015 ])) -1 (nil)
    (nil))

(note 491 490 492 63 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 492 491 493 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 68 [ __p.1009 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 493 492 494 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 500)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8100 [0x1fa4])
        (nil)))
;; End of basic block 63, registers live:
 (nil)

;; Start of basic block 64, registers live: (nil)
(code_label 494 493 495 64 107 "" [1 uses])

(note 495 494 496 64 [bb 64] NOTE_INSN_BASIC_BLOCK)

(note 496 495 497 64 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 497 496 498 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (mem/s:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                        (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                (plus:DI (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                            (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 498 497 499 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 506)) -1 (nil)
    (nil))
;; End of basic block 64, registers live:
 (nil)

(barrier 499 498 500)

;; Start of basic block 65, registers live: (nil)
(code_label 500 499 501 65 106 "" [1 uses])

(note 501 500 502 65 [bb 65] NOTE_INSN_BASIC_BLOCK)

(note 502 501 503 65 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 503 502 504 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 68 [ __p.1009 ]) [5 S8 A64])
        (reg:DI 62 [ kmer.1018 ])) -1 (nil)
    (nil))

(jump_insn 504 503 505 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (label_ref 494)) -1 (nil)
    (nil))
;; End of basic block 65, registers live:
 (nil)

(barrier 505 504 506)

;; Start of basic block 66, registers live: (nil)
(code_label 506 505 507 66 105 "" [2 uses])

(note 507 506 508 66 [bb 66] NOTE_INSN_BASIC_BLOCK)

(note 508 507 509 66 ("./DNAsequence.cc") 56)

(insn 509 508 510 66 ./DNAsequence.cc:56 (set (mem/s:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -208 [0xffffffffffffff30])) [0 str+0 S1 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 510 509 511 66 ("./DNAsequence.cc") 57)

(insn 511 510 512 66 ./DNAsequence.cc:57 (set (reg/f:DI 162)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 512 511 513 66 ./DNAsequence.cc:57 (set (reg:DI 163)
        (mem/f:DI (plus:DI (reg/f:DI 162)
                (const_int 16 [0x10])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 513 512 514 66 ./DNAsequence.cc:57 (parallel [
            (set (reg:DI 164)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 514 513 515 66 ./DNAsequence.cc:57 (set (reg:DI 4 si)
        (reg:DI 163)) -1 (nil)
    (nil))

(insn 515 514 516 66 ./DNAsequence.cc:57 (set (reg:DI 5 di)
        (reg:DI 164)) -1 (nil)
    (nil))

(call_insn 516 515 517 66 ./DNAsequence.cc:57 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 517 516 518 66 ("./DNAsequence.cc") 58)

(insn 518 517 519 66 ./DNAsequence.cc:58 (set (reg:DI 165)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 519 518 520 66 ./DNAsequence.cc:58 (set (reg:DI 166)
        (mem/f:DI (plus:DI (reg:DI 71 [ pretmp.977 ])
                (reg:DI 165)) [5 S8 A64])) -1 (nil)
    (nil))

(insn 520 519 521 66 ./DNAsequence.cc:58 (parallel [
            (set (reg:DI 167)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 521 520 522 66 ./DNAsequence.cc:58 (set (reg:DI 4 si)
        (reg:DI 166)) -1 (nil)
    (nil))

(insn 522 521 523 66 ./DNAsequence.cc:58 (set (reg:DI 5 di)
        (reg:DI 167)) -1 (nil)
    (nil))

(call_insn 523 522 524 66 ./DNAsequence.cc:58 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 524 523 525 66 ("./DNAsequence.cc") 59)

(insn 525 524 526 66 ./DNAsequence.cc:59 (parallel [
            (set (reg:DI 168)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 526 525 527 66 ./DNAsequence.cc:59 (set (reg:DI 5 di)
        (reg:DI 168)) -1 (nil)
    (nil))

(call_insn 527 526 528 66 ./DNAsequence.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 528 527 529 66 ./DNAsequence.cc:59 (set (reg/f:DI 169)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 169)
        (nil)))

(insn 529 528 530 66 ./DNAsequence.cc:59 (set (reg:DI 61 [ kmer.1019 ])
        (reg/f:DI 169)) -1 (nil)
    (nil))

(insn 530 529 531 66 ./DNAsequence.cc:59 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -40 [0xffffffffffffffd8])) [5 kmer+0 S8 A64])
        (reg:DI 61 [ kmer.1019 ])) -1 (nil)
    (nil))

(insn 531 530 532 66 ./DNAsequence.cc:59 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 61 [ kmer.1019 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 532 531 534 66 ./DNAsequence.cc:59 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 205)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 66, registers live:
 (nil)

;; Start of basic block 67, registers live: (nil)
(note 534 532 535 67 [bb 67] NOTE_INSN_BASIC_BLOCK)

(note 535 534 536 67 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 604)

(insn 536 535 537 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:DI 64 [ temp.1016 ])
        (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 537 536 538 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 64 [ temp.1016 ])
            (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                    (const_int 16 [0x10])) [4 <variable>.D.34734._M_impl._M_end_of_storage+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 538 537 540 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:604 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 550)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9087 [0x237f])
        (nil)))
;; End of basic block 67, registers live:
 (nil)

;; Start of basic block 68, registers live: (nil)
(note 540 538 541 68 [bb 68] NOTE_INSN_BASIC_BLOCK)

(note 541 540 542 68 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 610)

(insn 542 541 543 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 82 [ D.39442 ])
        (reg:DI 64 [ temp.1016 ])) -1 (nil)
    (nil))

(insn 543 542 544 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (parallel [
            (set (reg:DI 170)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -40 [0xffffffffffffffd8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 544 543 545 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 1 dx)
        (reg:DI 170)) -1 (nil)
    (nil))

(insn 545 544 546 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 4 si)
        (reg:DI 82 [ D.39442 ])) -1 (nil)
    (nil))

(insn 546 545 547 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg:DI 5 di)
        (reg/v/f:DI 109 [ kmers ])) -1 (nil)
    (nil))

(call_insn 547 546 608 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (call (mem:QI (symbol_ref/i:DI ("_ZNSt6vectorIPcSaIS0_EE13_M_insert_auxEN9__gnu_cxx17__normal_iteratorIPS0_S2_EERKS0_") [flags 0x1] <function_decl 0x2acec22e5500 _M_insert_aux>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 68, registers live:
 (nil)

;; Start of basic block 69, registers live: (nil)
(note 608 547 548 69 [bb 69] NOTE_INSN_BASIC_BLOCK)

(jump_insn 548 608 549 69 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (pc)
        (label_ref 568)) -1 (nil)
    (nil))
;; End of basic block 69, registers live:
 (nil)

(barrier 549 548 550)

;; Start of basic block 70, registers live: (nil)
(code_label 550 549 551 70 109 "" [1 uses])

(note 551 550 552 70 [bb 70] NOTE_INSN_BASIC_BLOCK)

(insn 552 551 553 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:610 (set (reg/v/f:DI 67 [ __p.1010 ])
        (reg:DI 64 [ temp.1016 ])) -1 (nil)
    (nil))

(note 553 552 554 70 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 554 553 555 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 67 [ __p.1010 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 555 554 556 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 562)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8100 [0x1fa4])
        (nil)))
;; End of basic block 70, registers live:
 (nil)

;; Start of basic block 71, registers live: (nil)
(code_label 556 555 557 71 113 "" [1 uses])

(note 557 556 558 71 [bb 71] NOTE_INSN_BASIC_BLOCK)

(note 558 557 559 71 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h") 607)

(insn 559 558 560 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (parallel [
            (set (mem/s:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                        (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                (plus:DI (mem/s/f:DI (plus:DI (reg/v/f:DI 109 [ kmers ])
                            (const_int 8 [0x8])) [4 <variable>.D.34734._M_impl._M_finish+0 S8 A64])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 560 559 561 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_vector.h:607 (set (pc)
        (label_ref 568)) -1 (nil)
    (nil))
;; End of basic block 71, registers live:
 (nil)

(barrier 561 560 562)

;; Start of basic block 72, registers live: (nil)
(code_label 562 561 563 72 112 "" [1 uses])

(note 563 562 564 72 [bb 72] NOTE_INSN_BASIC_BLOCK)

(note 564 563 565 72 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h") 104)

(insn 565 564 566 72 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (mem/f:DI (reg/v/f:DI 67 [ __p.1010 ]) [5 S8 A64])
        (reg:DI 61 [ kmer.1019 ])) -1 (nil)
    (nil))

(jump_insn 566 565 567 72 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ext/new_allocator.h:104 (set (pc)
        (label_ref 556)) -1 (nil)
    (nil))
;; End of basic block 72, registers live:
 (nil)

(barrier 567 566 568)

;; Start of basic block 73, registers live: (nil)
(code_label 568 567 569 73 111 "" [2 uses])

(note 569 568 570 73 [bb 73] NOTE_INSN_BASIC_BLOCK)

(note 570 569 571 73 ("./DNAsequence.cc") 56)

(insn 571 570 572 73 ./DNAsequence.cc:56 (set (mem/s:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -208 [0xffffffffffffff30])) [0 str+0 S1 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 572 571 573 73 ("./DNAsequence.cc") 57)

(insn 573 572 574 73 ./DNAsequence.cc:57 (set (reg/f:DI 171)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -64 [0xffffffffffffffc0])) [4 nucs.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 574 573 575 73 ./DNAsequence.cc:57 (set (reg:DI 172)
        (mem/f:DI (plus:DI (reg/f:DI 171)
                (const_int 24 [0x18])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 575 574 576 73 ./DNAsequence.cc:57 (parallel [
            (set (reg:DI 173)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 576 575 577 73 ./DNAsequence.cc:57 (set (reg:DI 4 si)
        (reg:DI 172)) -1 (nil)
    (nil))

(insn 577 576 578 73 ./DNAsequence.cc:57 (set (reg:DI 5 di)
        (reg:DI 173)) -1 (nil)
    (nil))

(call_insn 578 577 579 73 ./DNAsequence.cc:57 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 579 578 580 73 ("./DNAsequence.cc") 58)

(insn 580 579 581 73 ./DNAsequence.cc:58 (set (reg:DI 174)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -96 [0xffffffffffffffa0])) [4 km1mers.D.34734._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 581 580 582 73 ./DNAsequence.cc:58 (set (reg:DI 175)
        (mem/f:DI (plus:DI (reg:DI 71 [ pretmp.977 ])
                (reg:DI 174)) [5 S8 A64])) -1 (nil)
    (nil))

(insn 582 581 583 73 ./DNAsequence.cc:58 (parallel [
            (set (reg:DI 176)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 583 582 584 73 ./DNAsequence.cc:58 (set (reg:DI 4 si)
        (reg:DI 175)) -1 (nil)
    (nil))

(insn 584 583 585 73 ./DNAsequence.cc:58 (set (reg:DI 5 di)
        (reg:DI 176)) -1 (nil)
    (nil))

(call_insn 585 584 586 73 ./DNAsequence.cc:58 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strcat") [flags 0x41] <function_decl 0x2acebfb0ed00 strcat>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 586 585 587 73 ("./DNAsequence.cc") 59)

(insn 587 586 588 73 ./DNAsequence.cc:59 (parallel [
            (set (reg:DI 177)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -208 [0xffffffffffffff30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 588 587 589 73 ./DNAsequence.cc:59 (set (reg:DI 5 di)
        (reg:DI 177)) -1 (nil)
    (nil))

(call_insn 589 588 590 73 ./DNAsequence.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2acebfb11c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 590 589 591 73 ./DNAsequence.cc:59 (set (reg/f:DI 178)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 178)
        (nil)))

(insn 591 590 592 73 ./DNAsequence.cc:59 (set (reg:DI 63 [ kmer.1017 ])
        (reg/f:DI 178)) -1 (nil)
    (nil))

(insn 592 591 593 73 ./DNAsequence.cc:59 (set (mem/f/c/i:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -40 [0xffffffffffffffd8])) [5 kmer+0 S8 A64])
        (reg:DI 63 [ kmer.1017 ])) -1 (nil)
    (nil))

(insn 593 592 594 73 ./DNAsequence.cc:59 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 63 [ kmer.1017 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 594 593 597 73 ./DNAsequence.cc:59 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 205)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 73, registers live:
 (nil)

;; Start of basic block 74, registers live: (nil)
(note 597 594 595 74 [bb 74] NOTE_INSN_BASIC_BLOCK)

(jump_insn 595 597 596 74 ./DNAsequence.cc:59 (set (pc)
        (label_ref 224)) -1 (nil)
    (nil))
;; End of basic block 74, registers live:
 (nil)

(barrier 596 595 598)

(note 598 596 599 NOTE_INSN_FUNCTION_END)

(note 599 598 603 ("./DNAsequence.cc") 67)

;; Start of basic block 75, registers live: (nil)
(code_label 603 599 610 75 114 "" [4 uses])

(note 610 603 0 75 [bb 75] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 75, registers live:
 (nil)

