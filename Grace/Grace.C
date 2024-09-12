/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Grace.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Grace, 0);
    addToRunTimeSelectionTable(dragModel, Grace, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Grace::Grace
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Grace::~Grace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Grace::CdRe() const
{
    volScalarField Re(pair_.Re());
    volScalarField Eo(pair_.Eo());
    volScalarField dd(pair_.dispersed().d());

    volScalarField mud(pair_.dispersed().mu());
    volScalarField muc(pair_.continuous().mu());

    volScalarField rhoc(pair_.continuous().rho());
    volScalarField rhod(pair_.dispersed().rho());

    volScalarField Mo(pair_.Mo());
    dimensionedScalar muref
    (
        "muref",
        dimensionSet(1, -1, -1, 0, 0),
        scalar(0.0009)
    );

    volScalarField H
    (
        4.0/3.0*Eo*pow(Mo, -0.149)*pow(muc/muref, -0.14)
    );

    //Info<<"t"<<H;

    volScalarField J
    (
        neg(H - 59.3)*0.94*pow(H, 0.751)
      + pos(H - 59.3)*3.42*pow(H, 0.441)
    );

    //Info<<"Debug "<<dd;

    dimensionedScalar grav
    (
        "grav",
        dimensionSet(0, 1, -2, 0, 0),
        scalar(9.81)
    );

    //Info <<"Debug" << 1.0/dd;

    volScalarField Ut
    (
        muc/rhoc/dd*pow(Mo, -0.149)*(J - 0.857)
    );

    //Info <<"Debug" << nl;

    return 
        4.0/3.0*(rhoc - rhod)/rhoc
        *grav*dd/sqr(Ut)
        *max(Re, residualRe_);
}


// ************************************************************************* //
