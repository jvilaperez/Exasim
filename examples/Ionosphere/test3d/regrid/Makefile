OBJS := $(addsuffix .o, $(basename $(shell cat Srcfiles)))

$(EXEC): $(OBJS)
	$(F90) -o $@ $^ $(LDFLAGS)

%.o: %.F90
	$(F90) -c $(FFLAGS) $<

%.o: %.F
	$(F90) -c $(FFLAGS) $<

include Depends

all: $(EXEC)

clean:
	rm $(OBJS) *.mod Depends Srcfiles
