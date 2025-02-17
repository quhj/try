#include "readGravitationalAcceleration.H"

word continuousPhaseName
(
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    ).lookup("continuousPhaseName")
);

Info<< "Reading field U\n" << endl;
volVectorField Uc
(
    IOobject
    (
        IOobject::groupName("U", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field pOverRho\n" << endl;
volScalarField pOverRho
(
    IOobject
    (
        "pOverRho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating continuous-phase face flux field phic\n"
    << endl;

surfaceScalarField phic
(
    IOobject
    (
        IOobject::groupName("phi", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(pOverRho, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(pOverRho.name());

Info<< "Creating turbulence model\n" << endl;

singlePhaseTransportModel continuousPhaseTransport(Uc, phic);

dimensionedScalar rhocValue
(
    IOobject::groupName("rho", continuousPhaseName),
    dimDensity,
    continuousPhaseTransport.lookup
    (
        IOobject::groupName("rho", continuousPhaseName)
    )
);

volScalarField rhoc
(
    IOobject
    (
        rhocValue.name(),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhocValue
);

volScalarField rhoUnit 
(
    IOobject
    (
        "rhoUnit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar
	(
		"rhoUnit",
		dimensionSet(1, -3, 0, 0, 0),
		1.0
	)
);

volScalarField muc
(
    IOobject
    (
        IOobject::groupName("mu", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoc*continuousPhaseTransport.nu()
);

Info << "Creating field alphac\n" << endl;
// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alphac
(
    IOobject
    (
        IOobject::groupName("alpha", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

word kinematicCloudName("kinematicCloud");
args.optionReadIfPresent("cloudName", kinematicCloudName);

Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
basicKinematicTypeCloud kinematicCloud
(
    kinematicCloudName,
    rhoc,
    Uc,
    muc,
    g
);

// Particle fraction upper limit
scalar alphacMin
(
    1.0
  - readScalar
    (
        kinematicCloud.particleProperties().subDict("constantProperties")
       .lookup("alphaMax")
    )
);

Switch calculatePU
(
	kinematicCloud.particleProperties()
	.subDict("solution").lookup("calculatePU")
);

// Update alphac from the particle locations
alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));
surfaceScalarField alphaPhic("alphaPhic", alphacf*phic);

autoPtr<PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel>>
continuousPhaseTurbulence
(
    PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
    (
        alphac,
        Uc,
        alphaPhic,
        phic,
        continuousPhaseTransport
    )
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
	pOverRho*rhoUnit
);

uniformDimensionedScalarField hRef
(
	IOobject
	(
		"hRef",
		runTime.constant(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	dimensionedScalar
	(
		"hRef", 
		dimLength, 
		0.0
	)
);

dimensionedScalar ghRef
(
	mag(g.value()) > SMALL
  ? g & (cmptMag(g.value())/mag(g.value()))*hRef
  : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
);

Info<< "ghRef is " << ghRef.value() << nl;
volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

volVectorField gradP 
(
    IOobject
    (
		"gradP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	fvc::grad(p)
);
