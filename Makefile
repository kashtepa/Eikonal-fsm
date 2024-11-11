all: fsm

fsm:
	julia ./test/calc_fsm.jl luneburg
	julia ./test/calc_fsm.jl maxwell
	julia ./test/calc_fsm_eaton.jl

clean:
	-rm -vf img/*.pdf
	-rm -vf img/*.png

