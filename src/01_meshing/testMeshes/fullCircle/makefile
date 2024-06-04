all: clean domain

domain:
	gfortran -o buildDomain -std=f2008 -ffree-form buildDomain.f
	chmod +x buildDomain
	./buildDomain

clean:
	rm -f domain.poly
	rm -f buildDomain
