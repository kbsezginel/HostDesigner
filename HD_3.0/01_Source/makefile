FC     = gfortran
FFLAGS =  -ffixed-form -ffixed-line-length-72 -fcray-pointer -fbounds-check -ff2c 
SRC1 = alignz.f alignzpass.f angleval.f assignclass.f bondsplit.f\
       bumpcheckc.f bumpcheckd.f bumpcheckover.f center.f connlist.f\
       connlistg.f countlib.f drive.f driver.f equivrotor.f\
       fetchrecord.f getbondl.f getbondo.f getrmsd.f gettime.f\
       guestsym.f hitcheck.f hostdesigner.f inertia.f jacobi.f\
       linker.f makec.f maked.f makeover.f makeview.f makeview2.f\
       metalshape.f overlay.f quatfit.f readconstants.f readcontrol.f\
       readhosta.f readhostb.f readlink.f setbump.f setdihedral.f\
       sortarrays.f sortinteger.f sortreal.f summary.f superimpose.f\
       torsionval.f vieworient.f writedrives.f

OBJ1 = alignz.o alignzpass.o angleval.o assignclass.o bondsplit.o\
       bumpcheckc.o bumpcheckd.o bumpcheckover.o center.o connlist.o\
       connlistg.o countlib.o drive.o driver.o equivrotor.o\
       fetchrecord.o getbondl.o getbondo.o getrmsd.o gettime.o\
       guestsym.o hitcheck.o hostdesigner.o inertia.o jacobi.o\
       linker.o makec.o maked.o makeover.o makeview.o makeview2.o\
       metalshape.o overlay.o quatfit.o readconstants.o readcontrol.o\
       readhosta.o readhostb.o readlink.o setbump.o setdihedral.o\
       sortarrays.o sortinteger.o sortreal.o summary.o superimpose.o\
       torsionval.o vieworient.o writedrives.o

OBJS = $(OBJ1) 

hd3.0: $(OBJS)
	$(FC) $(FFLAGS) $(FOPT) $(LIBS) $(OBJS) -o hd3.0

.f.o:
	$(FC) $(FFLAGS) $(FOPT) -c $*.f

clean:
	rm -f *.o hd3.0
