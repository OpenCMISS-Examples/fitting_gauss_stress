#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to fit a stress field from the stress at Gauss points using OpenCMISS calls in python.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): Chris Bradley
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> Main script
# Add Python bindings directory to PATH
import sys, os

# Intialise OpenCMISS
from opencmiss.opencmiss import OpenCMISS_Python as oc

# Set problem parameters

# Cantilever dimensions
width = 10.0
length = 30.0
height = 10.0

force = -0.3

mooneyRivlin1 = 2.0
mooneyRivlin2 = 2.0
density = 9.0e-4

#pInit = -8.0 # Initial hydrostatic pressure
pInit = 0.0 # Initial hydrostatic pressure

pRef = pInit # Reference hydrostatic pressure

numberOfGaussXi = 3

numberOfLoadIncrements = 1

# Should not need to change anything below here

contextUserNumber = 1
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
decomposerUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
elasticityEquationsSetUserNumber = 1
elasticityEquationsSetFieldUserNumber = 3
elasticityDependentFieldUserNumber = 4
elasticityMaterialsFieldUserNumber = 5
elasticityStressFieldUserNumber = 6
elasticityProblemUserNumber = 1
fittingEquationsSetUserNumber = 2
fittingEquationsSetFieldUserNumber = 9
fittingDependentFieldUserNumber = 10
fittingMaterialsFieldUserNumber = 11
fittingProblemUserNumber = 2

numberOfGauss = pow(numberOfGaussXi,3)

if len(sys.argv) > 1:
    if len(sys.argv) > 6:
        sys.exit('Error: too many arguments- currently only accepting 5 options: numberXElements numberYElements numberZElements tau kappa')
    if len(sys.argv) >= 4:
        numberOfGlobalXElements = int(sys.argv[1])
        numberOfGlobalYElements = int(sys.argv[2])
        numberOfGlobalZElements = int(sys.argv[3])
    else:
        numberOfGlobalXElements = 1
        numberOfGlobalYElements = 1
        numberOfGlobalZElements = 3
    if len(sys.argv) == 6:
        tau = float(sys.argv[4])
        kappa = float(sys.argv[5])
    else:
        tau = 0.01
        kappa = 0.0005
else:
    numberOfGlobalXElements = 1
    numberOfGlobalYElements = 1
    numberOfGlobalZElements = 3
    tau = 0.001
    kappa = 0.00005

numberOfNodesXi = 3
numberOfXNodes = numberOfGlobalXElements*(numberOfNodesXi-1)+1
numberOfYNodes = numberOfGlobalYElements*(numberOfNodesXi-1)+1
numberOfZNodes = numberOfGlobalZElements*(numberOfNodesXi-1)+1
numberOfNodes = numberOfXNodes*numberOfYNodes*numberOfZNodes
    
context = oc.Context()
context.Create(contextUserNumber)

worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

# Get the number of computational nodes and this computational node number
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber,context)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = oc.Region()
region.CreateStart(regionUserNumber,worldRegion)
region.LabelSet("CantileverRegion")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define quadratic basis
quadraticBasis = oc.Basis()
quadraticBasis.CreateStart(basisUserNumber,context)
quadraticBasis.type = oc.BasisTypes.LAGRANGE_HERMITE_TP
quadraticBasis.numberOfXi = 3
quadraticBasis.interpolationXi = [oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*3
quadraticBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
quadraticBasis.CreateFinish()

# Define linear basis
linearBasis = oc.Basis()
linearBasis.CreateStart(pressureBasisUserNumber,context)
linearBasis.type = oc.BasisTypes.LAGRANGE_HERMITE_TP
linearBasis.numberOfXi = 3
linearBasis.interpolationXi = [oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
linearBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
linearBasis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = oc.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = oc.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [quadraticBasis,linearBasis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements]
# Finish the creation of a generated mesh in the region
mesh = oc.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = oc.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.CreateFinish()

# Decompose 
decomposer = oc.Decomposer()
decomposer.CreateStart(decomposerUserNumber,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

# Create a field for the geometry
geometricField = oc.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(oc.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(oc.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = oc.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(oc.FieldTypes.FIBRE)
fibreField.DecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(oc.FieldVariableTypes.U,"Fibre")
fibreField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,2)
fibreField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,2)
fibreField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,2)
fibreField.CreateFinish()

# Create the elasticity equations_set
elasticityEquationsSetField = oc.Field()
elasticityEquationsSet = oc.EquationsSet()
elasticityEquationsSetSpecification = [oc.EquationsSetClasses.ELASTICITY,
                                       oc.EquationsSetTypes.FINITE_ELASTICITY,
                                       oc.EquationsSetSubtypes.MOONEY_RIVLIN]
elasticityEquationsSet.CreateStart(elasticityEquationsSetUserNumber,region,fibreField,
                         elasticityEquationsSetSpecification,elasticityEquationsSetFieldUserNumber,
                         elasticityEquationsSetField)
elasticityEquationsSet.CreateFinish()

# Create the dependent field
elasticityDependentField = oc.Field()
elasticityEquationsSet.DependentCreateStart(elasticityDependentFieldUserNumber,elasticityDependentField)
elasticityDependentField.VariableLabelSet(oc.FieldVariableTypes.U,"ElasticityDependent")
elasticityDependentField.VariableLabelSet(oc.FieldVariableTypes.T,"ElasticityTraction")
# Set the pressure to be nodally based and use the second mesh component
elasticityDependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U,4,oc.FieldInterpolationTypes.NODE_BASED)
elasticityDependentField.ComponentInterpolationSet(oc.FieldVariableTypes.T,4,oc.FieldInterpolationTypes.NODE_BASED)
elasticityDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,4,2)
elasticityDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.T,4,2)
elasticityEquationsSet.DependentCreateFinish()

# Initialise elasticity dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,
    elasticityDependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1)
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,
    elasticityDependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2)
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3,
    elasticityDependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3)
oc.Field.ComponentValuesInitialiseDP(
    elasticityDependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,4,pInit)

# Create a field for the stress field. We will use this as the independent field for fitting to save a field copy. 
elasticityStressField = oc.Field()
elasticityStressField.CreateStart(elasticityStressFieldUserNumber,region)
elasticityStressField.TypeSet(oc.FieldTypes.GENERAL)
elasticityStressField.DecompositionSet(decomposition)
elasticityStressField.GeometricFieldSet(geometricField)
elasticityStressField.DependentTypeSet(oc.FieldDependentTypes.DEPENDENT)
elasticityStressField.NumberOfVariablesSet(2)
elasticityStressField.VariableTypesSet([oc.FieldVariableTypes.U,oc.FieldVariableTypes.V])
elasticityStressField.VariableLabelSet(oc.FieldVariableTypes.U,"GaussStress")
elasticityStressField.VariableLabelSet(oc.FieldVariableTypes.V,"GaussWeight")
elasticityStressField.NumberOfComponentsSet(oc.FieldVariableTypes.U,6)
elasticityStressField.NumberOfComponentsSet(oc.FieldVariableTypes.V,6)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,4,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,5,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,6,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.V,1,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.V,2,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.V,3,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.V,4,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.V,5,1)
elasticityStressField.ComponentMeshComponentSet(oc.FieldVariableTypes.V,6,1)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.U,1,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.U,2,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.U,3,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.U,4,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.U,5,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.U,6,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.V,1,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.V,2,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.V,3,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.V,4,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.V,5,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(oc.FieldVariableTypes.V,6,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.CreateFinish()

# Create the derived equations set stress fields
elasticityEquationsSet.DerivedCreateStart(elasticityStressFieldUserNumber,elasticityStressField)
elasticityEquationsSet.DerivedVariableSet(oc.EquationsSetDerivedTensorTypes.CAUCHY_STRESS,oc.FieldVariableTypes.U)
elasticityEquationsSet.DerivedCreateFinish()

# Create the material field
elasticityMaterialsField = oc.Field()
elasticityEquationsSet.MaterialsCreateStart(elasticityMaterialsFieldUserNumber,elasticityMaterialsField)
elasticityMaterialsField.VariableLabelSet(oc.FieldVariableTypes.U,"Material")
elasticityMaterialsField.VariableLabelSet(oc.FieldVariableTypes.V,"Density")
elasticityEquationsSet.MaterialsCreateFinish()

# Set materials parameters
elasticityMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES, \
                                                     1,mooneyRivlin1)
elasticityMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES, \
                                                     2,mooneyRivlin2)
elasticityMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES, \
                                                     1,density)

# Create elasticity equations
elasticityEquations = oc.Equations()
elasticityEquationsSet.EquationsCreateStart(elasticityEquations)
elasticityEquations.sparsityType = oc.EquationsSparsityTypes.SPARSE
elasticityEquations.outputType = oc.EquationsOutputTypes.NONE
elasticityEquationsSet.EquationsCreateFinish()

# Create the fitting equations_set
fittingEquationsSetField = oc.Field()
fittingEquationsSet = oc.EquationsSet()
fittingEquationsSetSpecification = [oc.EquationsSetClasses.FITTING,
                                    oc.EquationsSetTypes.GAUSS_FITTING_EQUATION,
                                    oc.EquationsSetSubtypes.GENERALISED_GAUSS_FITTING,
                                    oc.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
fittingEquationsSet.CreateStart(fittingEquationsSetUserNumber,region,geometricField,
                                fittingEquationsSetSpecification,fittingEquationsSetFieldUserNumber,
                                fittingEquationsSetField)
fittingEquationsSet.CreateFinish()

# Create the fitting dependent field
fittingDependentField = oc.Field()
fittingEquationsSet.DependentCreateStart(fittingDependentFieldUserNumber,fittingDependentField)
fittingDependentField.VariableLabelSet(oc.FieldVariableTypes.U,"FittedStress")
# Set the number of components to 2
fittingDependentField.NumberOfComponentsSet(oc.FieldVariableTypes.U,6)
# Set the field variables to be triquadratic Lagrange
fittingDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
fittingDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
fittingDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
fittingDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,4,1)
fittingDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,5,1)
fittingDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,6,1)
# Finish creating the fitting dependent field
fittingEquationsSet.DependentCreateFinish()

# Create the fitting independent field. Use the previously created elasticity stress field
fittingEquationsSet.IndependentCreateStart(elasticityStressFieldUserNumber,elasticityStressField)
# Finish creating the fitting independent field
fittingEquationsSet.IndependentCreateFinish()

# Initialise Gauss point weight field to 1.0
elasticityStressField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,1,1.0)
elasticityStressField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,2,1.0)
elasticityStressField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,3,1.0)
elasticityStressField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,4,1.0)
elasticityStressField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,5,1.0)
elasticityStressField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,6,1.0)

# Create material field (Sobolev parameters)
fittingMaterialField = oc.Field()
fittingEquationsSet.MaterialsCreateStart(fittingMaterialsFieldUserNumber,fittingMaterialField)
fittingMaterialField.VariableLabelSet(oc.FieldVariableTypes.U,"SmoothingParameters")
fittingEquationsSet.MaterialsCreateFinish()
# Set kappa and tau - Sobolev smoothing parameters
fittingMaterialField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,tau)
fittingMaterialField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,kappa)

# Create the fitting equations
fittingEquations = oc.Equations()
fittingEquationsSet.EquationsCreateStart(fittingEquations)
# Set the fitting equations sparsity type
fittingEquations.sparsityType = oc.EquationsSparsityTypes.SPARSE
# Set the fitting equations output type to none
fittingEquations.outputType = oc.EquationsOutputTypes.NONE
# Finish creating the fitting equations
fittingEquationsSet.EquationsCreateFinish()

# Define the elasticity problem
elasticityProblem = oc.Problem()
elasticityProblemSpecification = [oc.ProblemClasses.ELASTICITY,
                                  oc.ProblemTypes.FINITE_ELASTICITY,
                                  oc.ProblemSubtypes.STATIC_FINITE_ELASTICITY]
elasticityProblem.CreateStart(elasticityProblemUserNumber,context,elasticityProblemSpecification)
elasticityProblem.CreateFinish()

# Create the elasticity problem control loop
elasticityProblem.ControlLoopCreateStart()
controlLoop = oc.ControlLoop()
elasticityProblem.ControlLoopGet([oc.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
elasticityProblem.ControlLoopCreateFinish()

# Create elasticity problem solvers
elasticityNonLinearSolver = oc.Solver()
elasticityLinearSolver = oc.Solver()
elasticityProblem.SolversCreateStart()
elasticityProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,elasticityNonLinearSolver)
elasticityNonLinearSolver.outputType = oc.SolverOutputTypes.MONITOR
#elasticityNonLinearSolver.outputType = oc.SolverOutputTypes.PROGRESS
#elasticityNonLinearSolver.outputType = oc.SolverOutputTypes.MATRIX
elasticityNonLinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.FD)
elasticityNonLinearSolver.NewtonAbsoluteToleranceSet(1e-14)
elasticityNonLinearSolver.NewtonSolutionToleranceSet(1e-14)
elasticityNonLinearSolver.NewtonRelativeToleranceSet(1e-14)
elasticityNonLinearSolver.NewtonLinearSolverGet(elasticityLinearSolver)
#elasticityLinearSolver.linearType = oc.LinearSolverTypes.DIRECT
elasticityProblem.SolversCreateFinish()

# Create elasticity solver equations and add elasticity equations set to solver equations
elasticitySolverEquations = oc.SolverEquations()
elasticityProblem.SolverEquationsCreateStart()
elasticityNonLinearSolver.SolverEquationsGet(elasticitySolverEquations)
elasticitySolverEquations.sparsityType = oc.SolverEquationsSparsityTypes.SPARSE
elasticityEquationsSetIndex = elasticitySolverEquations.EquationsSetAdd(elasticityEquationsSet)
elasticityProblem.SolverEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
elasticityBoundaryConditions = oc.BoundaryConditions()
elasticitySolverEquations.BoundaryConditionsCreateStart(elasticityBoundaryConditions)

for widthNodeIdx in range(1,numberOfXNodes+1):
    for heightNodeIdx in range(1,numberOfYNodes+1):
        # Set left hand build in nodes ot no displacement
        nodeIdx=widthNodeIdx+(heightNodeIdx-1)*numberOfXNodes
        elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,1,nodeIdx,1,
                                             oc.BoundaryConditionsTypes.FIXED,0.0)
        elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,1,nodeIdx,2,
                                             oc.BoundaryConditionsTypes.FIXED,0.0)
        elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,1,nodeIdx,3,
                                             oc.BoundaryConditionsTypes.FIXED,0.0)
        print("Wall node number = %d" % (nodeIdx))
    # Set downward force on right-hand edge
    nodeIdx=numberOfNodes-widthNodeIdx+1
    elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.T,1,1,nodeIdx,2,
                                         oc.BoundaryConditionsTypes.NEUMANN_POINT,force)
    print("Force node number = %d" % (nodeIdx))
# Set reference pressure
elasticityBoundaryConditions.AddNode(elasticityDependentField,oc.FieldVariableTypes.U,1,1,numberOfNodes,4,
                                     oc.BoundaryConditionsTypes.FIXED,pRef)

elasticitySolverEquations.BoundaryConditionsCreateFinish()

# Create fitting problem
fittingProblemSpecification = [oc.ProblemClasses.FITTING,
                               oc.ProblemTypes.FITTING,
                               oc.ProblemSubtypes.STATIC_LINEAR_FITTING]
fittingProblem = oc.Problem()
fittingProblem.CreateStart(fittingProblemUserNumber,context,fittingProblemSpecification)
fittingProblem.CreateFinish()

# Create control loops
fittingProblem.ControlLoopCreateStart()
fittingProblem.ControlLoopCreateFinish()

# Create problem solver
fittingSolver = oc.Solver()
fittingProblem.SolversCreateStart()
fittingProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,fittingSolver)
fittingSolver.outputType = oc.SolverOutputTypes.PROGRESS
fittingProblem.SolversCreateFinish()

# Create fitting solver equations and add fitting equations set to solver equations
fittingSolverEquations = oc.SolverEquations()
fittingProblem.SolverEquationsCreateStart()
# Get the solver equations
fittingSolver.SolverEquationsGet(fittingSolverEquations)
fittingSolverEquations.sparsityType = oc.SolverEquationsSparsityTypes.SPARSE
fittingEquationsSetIndex = fittingSolverEquations.EquationsSetAdd(fittingEquationsSet)
fittingProblem.SolverEquationsCreateFinish()

# Prescribe boundary conditions for the fitting problem
fittingBoundaryConditions = oc.BoundaryConditions()
fittingSolverEquations.BoundaryConditionsCreateStart(fittingBoundaryConditions)
fittingSolverEquations.BoundaryConditionsCreateFinish()

# Solve the elasticity problem
elasticityProblem.Solve()

# Calculate the stress field
elasticityEquationsSet.DerivedVariableCalculate(oc.EquationsSetDerivedTensorTypes.CAUCHY_STRESS)
                        
# Solve the fitting problem
fittingProblem.Solve()
                    
# Export results
fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("Cantilever","FORTRAN")
fields.ElementsExport("Cantilever","FORTRAN")
fields.Finalise()

