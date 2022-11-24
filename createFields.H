#include "createRDeltaT.H"

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
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Linearly independent tensors T0, ..., T9 made from
// different combinations of strain rate and rotation rate tensors.
volTensorField gradU
(
    IOobject
    (
        "gradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        //IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    fvc::grad(U)
);

volScalarField inv1GradU
(
    IOobject
    (
        "inv1GradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        //IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    tr(gradU)
);

volScalarField inv2GradU
(
    IOobject
    (
        "inv2GradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        //IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    tr(gradU) * tr(gradU) - tr(gradU & gradU)
);

volSymmTensorField strainRateTensor
(
    IOobject
    (
        "strainRateTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        //IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    symm(gradU)
);

volTensorField rotationRateTensor
(
    IOobject
    (
        "rotationRateTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
//        IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    skew(gradU)
);

volTensorField T0
(
    IOobject
    (
        "T0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tensor::I & strainRateTensor
);

volTensorField T1
(
    IOobject
    (
        "T1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (strainRateTensor & rotationRateTensor) - (rotationRateTensor & strainRateTensor)
);

volTensorField T2
(
    IOobject
    (
        "T2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (strainRateTensor & strainRateTensor) - (tensor::I * (tr(strainRateTensor & strainRateTensor) / 3))
);

volTensorField T3
(
    IOobject
    (
        "T3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rotationRateTensor & rotationRateTensor) - (tensor::I * (tr(rotationRateTensor & rotationRateTensor) / 3))
);

volTensorField T4
(
    IOobject
    (
        "T4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rotationRateTensor & strainRateTensor & strainRateTensor) -
		(strainRateTensor & strainRateTensor & rotationRateTensor)
);

volTensorField T5
(
    IOobject
    (
        "T5",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rotationRateTensor & rotationRateTensor & strainRateTensor) +
		(strainRateTensor & rotationRateTensor & rotationRateTensor) -
		(tensor::I * (2 * tr(strainRateTensor & rotationRateTensor & rotationRateTensor) / 3))
);

volTensorField T6
(
    IOobject
    (
        "T6",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rotationRateTensor & strainRateTensor & rotationRateTensor & rotationRateTensor) -
		(rotationRateTensor & rotationRateTensor & strainRateTensor & rotationRateTensor)
);

volTensorField T7
(
    IOobject
    (
        "T7",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (strainRateTensor & rotationRateTensor & strainRateTensor & strainRateTensor) -
		(strainRateTensor & strainRateTensor & rotationRateTensor & strainRateTensor)
);

volTensorField T8
(
    IOobject
    (
        "T8",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rotationRateTensor & rotationRateTensor & strainRateTensor & strainRateTensor) +
		(strainRateTensor & strainRateTensor & rotationRateTensor & rotationRateTensor) -
		(tensor::I * (2 * tr(strainRateTensor & strainRateTensor & rotationRateTensor & rotationRateTensor) / 3))
);

volTensorField T9
(
    IOobject
    (
        "T9",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rotationRateTensor & strainRateTensor & strainRateTensor & rotationRateTensor & rotationRateTensor) -
		(rotationRateTensor & rotationRateTensor & strainRateTensor & strainRateTensor & rotationRateTensor)
);

volScalarField I0
(
    IOobject
    (
        "I0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tr(strainRateTensor & strainRateTensor)
);

volScalarField I1
(
    IOobject
    (
        "I1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tr(rotationRateTensor & rotationRateTensor)
);

volScalarField I2
(
    IOobject
    (
        "I2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tr(strainRateTensor & strainRateTensor & strainRateTensor)
);

volScalarField I3
(
    IOobject
    (
        "I3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tr(rotationRateTensor & rotationRateTensor & strainRateTensor)
);

volScalarField I4
(
    IOobject
    (
        "I4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tr(rotationRateTensor & rotationRateTensor & strainRateTensor & strainRateTensor)
);

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

volScalarField magGradP
(
    IOobject
    (
        "magGradP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(gradP)
);

#include "createMRF.H"
#include "createFvOptions.H"
