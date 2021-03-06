CC	= g++
CFLAGS	= -Wall -lm
SDIR	= src
ODIR	= build
OUT	= out
TARGET	= bin/runner

# need to make a choice of ONE suffix
SRCEXT	= c

SRC	= $(wildcard $(SDIR)/*.$(SRCEXT))
OBJ	= $(patsubst $(SDIR)/%,$(ODIR)/%,$(SRC:.$(SRCEXT)=.o))
INC	= -I include

$(TARGET): $(OBJ)
	@mkdir -p bin
	@mkdir -p $(OUT)/data
	@echo "Compiling : $(CC) $^ -o $(TARGET)"; $(CC) $^ -o $(TARGET)

$(ODIR)/%.o: $(SDIR)/%.$(SRCEXT)
	@mkdir -p $(ODIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo "Cleaning : $(RM) -r $(ODIR) $(TARGET)"; $(RM) -r $(ODIR) $(TARGET); $(RM) bin/*;

# auxiliary compiles go here:

plotter1:
	gle -o "out/Z2_D4.pdf" -d pdf "out/plotter/Z2_D4.gle"
ZN:
	gle -o "out/Zn_D4.pdf" -d pdf "out/plotter/Zn_D4.gle"
plotter2:
	gle -o "out/compare_DIM.pdf" -d pdf "out/plotter/compare_DIM.gle"
plotter3:
	gle -o "out/U(1)_D4.pdf" -d pdf "out/plotter/U(1)_D4.gle"
plotter4:
	gle -o "out/tr_SU(N).pdf" -d pdf "out/plotter/tr_U.gle"
plotter5:
	gle -o "out/SU(N)_D4.pdf" -d pdf "out/plotter/SU(N)_D4.gle"

# test
tester: $(OBJ)
	$(CC) $(CFLAGS) $(INC) build/matrix.o spike/mat_test.c $(BIN) -o bin/tester

.PHONY: clean

