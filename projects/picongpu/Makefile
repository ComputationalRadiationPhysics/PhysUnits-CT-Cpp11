RM = rm

INCDIR = ../../
HDRDIR = $(INCDIR)/phys/units/
SRCDIR = ./

HEADERS = \
	io.hpp \
	io_output.hpp \
	io_output_eng.hpp \
	io_symbols.hpp \
	other_units.hpp \
	physical_constants.hpp \
	quantity.hpp \
	quantity_io.hpp \
	quantity_io_ampere.hpp \
	quantity_io_becquerel.hpp \
	quantity_io_candela.hpp \
	quantity_io_celsius.hpp \
	quantity_io_coulomb.hpp \
	quantity_io_dimensionless.hpp \
	quantity_io_engineering.hpp \
	quantity_io_farad.hpp \
	quantity_io_gray.hpp \
	quantity_io_henry.hpp \
	quantity_io_hertz.hpp \
	quantity_io_joule.hpp \
	quantity_io_kelvin.hpp \
	quantity_io_kilogram.hpp \
	quantity_io_lumen.hpp \
	quantity_io_lux.hpp \
	quantity_io_meter.hpp \
	quantity_io_mole.hpp \
	quantity_io_newton.hpp \
	quantity_io_ohm.hpp \
	quantity_io_pascal.hpp \
	quantity_io_radian.hpp \
	quantity_io_second.hpp \
	quantity_io_siemens.hpp \
	quantity_io_sievert.hpp \
	quantity_io_speed.hpp \
	quantity_io_steradian.hpp \
	quantity_io_symbols.hpp \
	quantity_io_tesla.hpp \
	quantity_io_volt.hpp \
	quantity_io_watt.hpp \
	quantity_io_weber.hpp \
	test_util.hpp

OBJS =

vpath %.hpp $(HDRDIR)
vpath %.cpp $(SRCDIR)

CC = g++
CXXFLAGS = -Wall -Wextra -Weffc++ -Wno-missing-braces -std=c++11 -DQUANTITY_USE_KELVIN -I$(INCDIR)

%.o: %.cpp $(HEADERS)

%.exe: %.o
	$(CC) -o $*.exe $<

all: main.exe run_tests

run_tests:
	./main.exe

clean:
	-$(RM) *.bak *.o

distclean: clean
	-$(RM) *.exe

#
# end of file
#
