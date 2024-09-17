from pathlib import Path

import json
import utils
import sherpa.astro.ui as shp
from gdpyc import GasMap
from sherpa.astro.instrument import RSPModelNoPHA
from sherpa.astro.xspec import read_xstable_model
from sherpa.models.parameter import Parameter


models_path = Path(".", "models")


def galactic_absorption(coords):
    nh = GasMap.nh(coords, nhmap='LAB')

    galabs = shp.xstbabs.galabs
    galabs.nH = nh.value / 1e22
    galabs.nH.freeze()

    return galabs


def background_only(bkgmodels):
    srcmodel = []
    parameters = []

    for i, bkg in enumerate(bkgmodels, 1):
        # diagonals rmf and arf
        rmf = bkg.parts[1].parts[0].rmf
        arf = bkg.parts[1].parts[0].arf

        norm = RSPModelNoPHA(arf, rmf, shp.xsconstant(f"bkgrenorm_{i}"))
        srcmodel.append(norm * bkg)

        lognorm = Parameter(
            modelname="src",
            name=f"log_bkgrenorm_{i}",
            val=0,
            min=-2,
            max=2,
            hard_min=-10,
            hard_max=10,
        )
        shp.link(f"bkgrenorm_{i}.factor", 10**lognorm)
        parameters.append(lognorm)

    return srcmodel, parameters


def _get_model_powerlaw(ids, fix_inst, use_tbabs_table=False, *args, **kwargs):
    if not use_tbabs_table:
        intabs = shp.xstbabs.intabs
    else:
        print("Using tbabs tabulated model")
        intabs = read_xstable_model("intabs", "/dataAGN/aetviita/WP6/xmm2athena_wp6/data/ttbabs.fits")
    po = shp.xspowerlaw.pl
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] *cfl(intabs * po))

    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=20, max=26, hard_min=20, hard_max=26
    )
    intabs.nh = 10**(lognH - 22)

    #    lognorm = Parameter(
    #        modelname="src", name="lognorm", val=0, min=-8, max=3, hard_min=-20, hard_max=20
    #    )
    #    po.norm = 10**lognorm

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = po.PhoIndex

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    param3.min = 0
    param3.max = 6


    parameters = [param1, param2, param3]

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters

def _get_model_powerlaw2(ids, fix_inst):
    intabs = shp.xstbabs.intabs
    po = shp.xspowerlaw.pl
    cfl = shp.xscflux.cflux
    const = shp.xsconstant.constant
    srcmodel = []

    inner_model = cfl(intabs * po)
    srcmodel.append(inner_model)
    if len(ids) > 1:
        srcmodel.append(const*inner_model)


    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=20, max=26, hard_min=20, hard_max=26
    )
    intabs.nh = 10**(lognH - 22)

    #    lognorm = Parameter(
    #        modelname="src", name="lognorm", val=0, min=-8, max=3, hard_min=-20, hard_max=20
    #    )
    #    po.norm = 10**lognorm

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = po.PhoIndex

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    param3.min = 0
    param3.max = 6


    parameters = [param1, param2, param3]


    #for i, id in enumerate(ids):
    #    if id == 1:
            #const[i].factor.val = 1.0
    #        const[i].factor.freeze()
    #    else:
    #        const[i].factor.min = 0.0
    #        const[i].factor.max = 2.0
    #        if fix_inst == 1:
    #            const[i].factor.freeze()
    #        else:
    #            parameters += [const[i].factor]

    if len(ids) > 1:
        const.factor.min = 0.0
        const.factor.max = 5.0
        if fix_inst == 1:
            const.factor.freeze()
        else:
            parameters += [const.factor]

    return srcmodel, parameters

def _get_model_plbestfit(ids, fix_inst, output_path):
    intabs = shp.xstbabs.intabs
    po = shp.xspowerlaw.pl
    cfl = shp.xscflux.cflux
    const = shp.xsconstant.constant
    srcmodel = []

    inner_model = cfl(intabs * po)
    srcmodel.append(inner_model)
    if len(ids) > 1:
        srcmodel.append(const*inner_model)


    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=20, max=26, hard_min=18, hard_max=28
    )
    intabs.nh = 10**(lognH - 22)

    #    lognorm = Parameter(
    #        modelname="src", name="lognorm", val=0, min=-8, max=3, hard_min=-20, hard_max=20
    #    )
    #    po.norm = 10**lognorm

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = po.PhoIndex

    info_path = output_path / "info"
    json_path = info_path.joinpath("results.json")
    with json_path.open("r") as rf:
        result = json.load(rf)

    bstpar = result['posterior']['mean']
    bstdev = result['posterior']['stdev']


    param1.val = bstpar[0]
    param1.min = bstpar[0]-bstdev[0]
    param1.max = bstpar[0]+bstdev[0]

    param2.val = bstpar[1]
    if bstpar[1]-bstdev[1] < 18:
        param2.min = 18
    else:
        param2.min = bstpar[1]-bstdev[1]

    if bstpar[1]+bstdev[1] > 28:
        param2.max = 28
    else:
        param2.max = bstpar[1]+bstdev[1]

    param3.val = bstpar[2]
    param3.min = bstpar[2]-bstdev[2]
    param3.max = bstpar[2]+bstdev[2]


    parameters = [param1, param2, param3]


    #for i, id in enumerate(ids):
    #    if id == 1:
            #const[i].factor.val = 1.0
    #        const[i].factor.freeze()
    #    else:
    #        const[i].factor.min = 0.0
    #        const[i].factor.max = 2.0
    #        if fix_inst == 1:
    #            const[i].factor.freeze()
    #        else:
    #            parameters += [const[i].factor]

    if len(ids) > 1:
        const.factor.val = bstpar[3]
        if bstpar[3]-bstdev[3] < 0:
            const.factor.min = 0
        else:
            const.factor.min = bstpar[3]-bstdev[3]

        const.factor.max = bstpar[3]+bstdev[3]
        if fix_inst == 1:
            const.factor.freeze()
        else:
            parameters += [const.factor]

    return srcmodel, parameters

def _get_model_powerlaw3(ids, fix_inst):
    intabs = shp.xstbnew.intabs
    po = shp.xspowerlaw.pl
    cfl = shp.xscflux.cflux
    const = shp.xsconstant.constant
    srcmodel = []

    inner_model = cfl(intabs * po)
    srcmodel.append(inner_model)
    if len(ids) > 1:
        srcmodel.append(const*inner_model)


    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=20, max=26, hard_min=20, hard_max=26
    )
    intabs.nh = 10**(lognH - 22)

    #    lognorm = Parameter(
    #        modelname="src", name="lognorm", val=0, min=-8, max=3, hard_min=-20, hard_max=20
    #    )
    #    po.norm = 10**lognorm

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = po.PhoIndex

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    param3.min = 0
    param3.max = 6


    parameters = [param1, param2, param3]


    #for i, id in enumerate(ids):
    #    if id == 1:
            #const[i].factor.val = 1.0
    #        const[i].factor.freeze()
    #    else:
    #        const[i].factor.min = 0.0
    #        const[i].factor.max = 2.0
    #        if fix_inst == 1:
    #            const[i].factor.freeze()
    #        else:
    #            parameters += [const[i].factor]

    if len(ids) > 1:
        const.factor.min = 0.0
        const.factor.max = 5.0
        if fix_inst == 1:
            const.factor.freeze()
        else:
            parameters += [const.factor]

    return srcmodel, parameters


def _get_model_blackbody(ids, fix_inst, *args, **kwargs):
    #intabs = shp.xstbabs.intabs
    intabs = shp.xsphabs.phabs
    bb = shp.xsbbody.bb
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * cfl(intabs * bb))

    lognH = Parameter(
            modelname="src", name="logNH", val=22, min=18, max=26, hard_min=18, hard_max=26
        )
    intabs.nh = 10**(lognH - 22)

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = bb.kT

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    param3.min = 0.01
    param3.max = 10

    parameters = [param1, param2, param3]


    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters


def _get_model_apec_single(ids, fix_inst, *args, **kwargs):
    intabs = shp.xstbabs.intabs
    apc = shp.xsapec.apc
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * cfl(intabs * apc))

    lognH = Parameter(
            #modelname="src", name="logNH", val=22, min=18, max=26, hard_min=18, hard_max=26
            modelname="src", name="logNH", val=22, min=18, max=24, hard_min=18, hard_max=24
        )
    intabs.nh = 10**(lognH - 22)

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = apc.kT

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    #param3.min = 0.01
    #param3.max = 20
    #param3.min = 0.10
    #param3.max = 1.70
    param3.min = 0.01
    param3.max = 17.0

    parameters = [param1, param2, param3]

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters


def _get_model_zpowlaw(ids, fix_inst, redshift=-1, *args, **kwargs):
    intabs = shp.xsztbabs.intabs
    po = shp.xszpowerlw.po
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * cfl(intabs * po))

    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=20, max=26, hard_min=20, hard_max=26
    )
    intabs.nh = 10 ** (lognH - 22)

    cfl.Emin = 0.2
    cfl.Emax = 12.0
    cfl.lg10Flux.min = -17
    cfl.lg10Flux.max = -7

    po.PhoIndex.min = 0
    po.PhoIndex.max = 6

    parameters = [cfl.lg10Flux, lognH, po.PhoIndex]

    # 20230131 AV: handle redshift information
    #  z <  0  redshift is a free parameter
    #  z >= 0  redshift is fixed
    if redshift < 0:
        raise ValueError("Negative redshift not expected. Check the priors.py file before trying this further. Regards, AV 20230131")
        po.redshift.thaw()
        intabs.redshift = po.redshift
        parameters += [po.redshift]
    else:
        intabs.redshift.val = redshift
        intabs.redshift.freeze()
        po.redshift.val = redshift
        po.redshift.freeze()
    print("Redshift set to", redshift)

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters

def _get_model_double_zpowlaw(ids, fix_inst, redshift=-1, *args, **kwargs):
    intabs = shp.xsztbabs.intabs
    po1 = shp.xszpowerlw.po1
    po2 = shp.xszpowerlw.po2
    cfl1 = shp.xscflux.cflux1
    cfl2 = shp.xscflux.cflux2
    #cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * (cfl1(intabs * po1) + cfl2(po2)))
        #srcmodel.append(const[i] * (cfl(intabs * (po1 + po2))))
	#srcmodel.append(const[i] * cfl(intabs * po1) + const[i] * po2))

    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=20, max=26, hard_min=20, hard_max=26
    )
    intabs.nh = 10 ** (lognH - 22)

    cfl1.Emin = 0.2
    cfl1.Emax = 12.0
    cfl1.lg10Flux.min = -17
    cfl1.lg10Flux.max = -7

    cfl2.Emin = 0.2
    cfl2.Emax = 12.0
    cfl2.lg10Flux.min = -17
    cfl2.lg10Flux.max = -7

    #cfl.Emin = 0.2
    #cfl.Emax = 12.0
    #cfl.lg10Flux.min = -17
    #cfl.lg10Flux.max = -7

    po1.PhoIndex.min = 0
    po1.PhoIndex.max = 6
    po2.PhoIndex.min = 0
    po2.PhoIndex.max = 6

    parameters = [cfl1.lg10Flux, lognH, po1.PhoIndex, cfl2.lg10Flux, po2.PhoIndex]
    #parameters = [cfl.lg10Flux, lognH, po1.PhoIndex, po2.PhoIndex]
   
    if redshift < 0:
        raise ValueError("Negative redshift not expected. Check the priors.py file before trying this further.")
        po1.redshift.thaw()
        intabs.redshift = po1.redshift
        po2.redshift = po1.redshift
        parameters += [po1.redshift]
    else:
        intabs.redshift.val = redshift
        intabs.redshift.freeze()
        po1.redshift.val = redshift
        po1.redshift.freeze()
        po2.redshift.val = redshift
        po2.redshift.freeze()

    print("Redshift set to", redshift)

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters


#def _get_model_zpowlaw(redshift):
#    intabs = shp.xsztbabs.intabs
#    po = shp.xszpowerlw.po
#    srcmodel = intabs * po

#    lognH = Parameter(
#        modelname="src", name="logNH", val=22, min=20, max=26, hard_min=20, hard_max=26
#    )
#    intabs.nh = 10**(lognH - 22)

#    lognorm = Parameter(
#        modelname="src", name="lognorm", val=0, min=-8, max=3, hard_min=-20, hard_max=20
#    )
#    po.norm = 10**lognorm

#    intabs.redshift = po.redshift

#    parameters = [
#        lognH,
#        po.PhoIndex,
#        lognorm,
#    ]

#    if redshift < 0:
#        po.redshift.thaw()
#        parameters += [po.redshift]
#    else:
#        po.redshift.val = redshift
#        po.redshift.freeze()

#    return srcmodel, parameters


def _get_model_torus(redshift):
    torus = read_xstable_model(
        "torus", models_path.joinpath("uxclumpy-cutoff.fits").as_posix()
    )
    scattering = read_xstable_model(
        "scattering", models_path.joinpath("uxclumpy-cutoff-omni.fits").as_posix()
    )
    srcmodel = torus + scattering

    lognH = Parameter(modelname="src", name="logNH", val=22, min=20, max=26, hard_min=20, hard_max=26)
    torus.nh = 10**(lognH - 22)
    scattering.nh = torus.nh

    # the limits correspond to fluxes between Sco X-1 and CDFS7Ms faintest fluxes
    torus_lognorm = Parameter(
        modelname="src", name="torus_lognorm", val=0, min=-8, max=3, hard_min=-20, hard_max=20
    )
    torus.norm = 10**torus_lognorm

    scattering_lognorm = Parameter(
        'src', 'scattering_lognorm', val=-2, min=-7, max=-1, hard_min=-7, hard_max=-1
    )
    scattering.norm = 10**(torus_lognorm + scattering_lognorm)
    #shp.set_par("scattering.norm", val=0, frozen=True)

    # Edge-on torus
    #shp.set_par("torus.theta_inc", val=45, frozen=False)
    #torus.theta_inc.val = 85

    scattering.phoindex = torus.phoindex
    scattering.ecut = torus.ecut
    scattering.theta_inc = torus.theta_inc
    scattering.torsigma = torus.torsigma
    scattering.ctkcover = torus.ctkcover
    scattering.redshift = torus.redshift
    scattering.redshift.max = 10

    parameters = [
        lognH,
        torus.phoindex,
        torus.theta_inc,
        torus_lognorm,
        scattering_lognorm
    ]

    if redshift < 0:
        torus.redshift.thaw()
        parameters += [torus.redshift]
    else:
        torus.redshift.val = redshift
        torus.redshift.freeze()

    torus.redshift.min = 0
    torus.redshift.max = 10

    return srcmodel, parameters


def _get_model_bremss(ids, fix_inst, *args, **kwargs):
    intabs = shp.xstbabs.intabs
    br = shp.xsbremss.mdl
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * cfl(intabs * br))

    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=18, max=24, hard_min=18, hard_max=24
    )
    intabs.nh = 10**(lognH - 22)

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = br.kT

    param1.min = -17
    param1.max = -7

    param3.min = 0.0001
    param3.max = 200.00

    parameters = [param1, param2, param3]

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters


def _get_model_apec_apec(ids, fix_inst, *args, **kwargs):
    intabs = shp.xstbabs.intabs
    apec1 = shp.xsapec("apec1")
    apec2 = shp.xsapec("apec2")
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * cfl(intabs * (apec1 + apec2)))

    lognH = Parameter(
            modelname="src", name="logNH", val=22, min=18, max=24, hard_min=18, hard_max=24
        )
    intabs.nh = 10**(lognH - 22)

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = apec1.kT
    param4 = apec2.kT

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    # The two APEC temperatures based on Ada Nebot email 2023-02-06
    param3.min = 0.10
    param3.max = 0.25

    param4.min = 0.30
    param4.max = 1.70

    parameters = [param1, param2, param3, param4]

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters

def _get_model_apec_apec_const(ids, fix_inst, *args, **kwargs):
    intabs = shp.xstbabs.intabs
    apec1 = shp.xsapec("apec1")
    apec2 = shp.xsapec("apec2")
    constant_apec = shp.xsconstant(f"constant_apec")
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * cfl(intabs * (apec1 + apec2 * constant_apec)))

    # Ada Nebot email 20230214: fix NH to 18-24
    lognH = Parameter(
            modelname="src", name="logNH", val=22, min=18, max=24, hard_min=18, hard_max=24
        )
    intabs.nh = 10**(lognH - 22)

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = apec1.kT
    param4 = apec2.kT
    param5 = constant_apec.factor

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    # The two APEC temperatures based on Ada Nebot email 2023-02-06
    param3.min = 0.10
    param3.max = 0.25

    param4.min = 0.30
    param4.max = 1.70

    param5.min = 0.
    param5.max = 100.

    parameters = [param1, param2, param3, param4, param5]

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters

def _get_model_powerlaw_blackbody(ids, fix_inst, *args, **kwargs):
    intabs = shp.xstbabs.intabs
    po = shp.xspowerlaw.pl
    bb = shp.xsbbody.bb
    cfl = shp.xscflux.cflux
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * cfl(intabs * (po + bb)))

    lognH = Parameter(
        modelname="src", name="logNH", val=22, min=18, max=24, hard_min=18, hard_max=24
    )
    intabs.nh = 10**(lognH - 22)

    #    lognorm = Parameter(
    #        modelname="src", name="lognorm", val=0, min=-8, max=3, hard_min=-20, hard_max=20
    #    )
    #    po.norm = 10**lognorm

    cfl.Emin = 0.2
    cfl.Emax = 12.0

    param1 = cfl.lg10Flux
    param2 = lognH
    param3 = po.PhoIndex
    param4 = bb.kT

    param1.min = -17
    param1.max = -7

    #param2.min = 0.01
    #param2.max = 1e4

    param3.min = 0
    param3.max = 6

    param4.min = 0.01
    param4.max = 10

    parameters = [param1, param2, param3, param4]

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters


def _get_model_powerlaw_blackbody_const(ids, fix_inst, *args, **kwargs):

    intabs = shp.xstbabs.intabs
    po = shp.xspowerlaw.pl
    bb = shp.xsbbody.bb
    cfl1 = shp.xscflux.cflux1
    cfl2 = shp.xscflux.cflux2
    const = []
    srcmodel = []

    for i, id in enumerate(ids):
        const.append(shp.xsconstant(f"constant_{id}"))
        srcmodel.append(const[i] * (cfl1(intabs * po) + cfl2(intabs * bb)))

    lognH = Parameter(modelname="src", name="logNH", val=22, min=18, max=24, hard_min=18, hard_max=24)
    intabs.nh = 10**(lognH - 22)

    cfl1.Emin = 0.2 ; cfl1.Emax = 12.0
    cfl2.Emin = 0.2 ; cfl2.Emax = 12.0

    param1 = cfl1.lg10Flux
    param2 = po.PhoIndex
    param3 = cfl2.lg10Flux
    param4 = bb.kT
    param5 = lognH

    param1.min = -17  ; param1.max = -7
    param2.min = 0    ; param2.max = 6

    param3.min = -17  ; param3.max = -7
    param4.min = 0.01 ; param4.max = 10

    parameters = [param1, param2, param3, param4, param5]

    if len(ids) == 1:
        const[0].factor.freeze()
    else:
        for i, id in enumerate(ids):
            detector = utils.get_detector(id)

            if detector == "EPN" or detector == "EPNA":
                const[i].factor.freeze()
            else:
                const[i].factor.min = 0.0
                const[i].factor.max = 5.0
                if fix_inst == 1:
                    const[i].factor.freeze()
                else:
                    parameters += [const[i].factor]

    return srcmodel, parameters


def get_source_model(model, *args, **kwargs):
    try:
        model = globals()[f"_get_model_{model}"]

    except KeyError:
        raise ValueError(f"Model '{model}' is not defined!!!")

    return model(*args, **kwargs)
