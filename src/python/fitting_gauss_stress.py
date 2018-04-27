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

# Intialise OpenCMISS-Iron
from opencmiss.iron import iron

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

coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
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
    
worldRegion = iron.Region()
iron.Context.WorldRegionGet(worldRegion)

# Get the number of computational nodes and this computational node number
computationEnvironment = iron.ComputationEnvironment()
iron.Context.ComputationEnvironmentGet(computationEnvironment)
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber,iron.Context)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,worldRegion)
region.LabelSet("CantileverRegion")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define quadratic basis
quadraticBasis = iron.Basis()
quadraticBasis.CreateStart(basisUserNumber,iron.Context)
quadraticBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
quadraticBasis.numberOfXi = 3
quadraticBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*3
quadraticBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
quadraticBasis.CreateFinish()

# Define linear basis
linearBasis = iron.Basis()
linearBasis.CreateStart(pressureBasisUserNumber,iron.Context)
linearBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
linearBasis.numberOfXi = 3
linearBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
linearBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
linearBasis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [quadraticBasis,linearBasis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements]
# Finish the creation of a generated mesh in the region
mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)
fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)
fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,2)
fibreField.CreateFinish()

# Create the elasticity equations_set
elasticityEquationsSetField = iron.Field()
elasticityEquationsSet = iron.EquationsSet()
elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                       iron.EquationsSetTypes.FINITE_ELASTICITY,
                                       iron.EquationsSetSubtypes.MOONEY_RIVLIN]
elasticityEquationsSet.CreateStart(elasticityEquationsSetUserNumber,region,fibreField,
                         elasticityEquationsSetSpecification,elasticityEquationsSetFieldUserNumber,
                         elasticityEquationsSetField)
elasticityEquationsSet.CreateFinish()

# Create the dependent field
elasticityDependentField = iron.Field()
elasticityEquationsSet.DependentCreateStart(elasticityDependentFieldUserNumber,elasticityDependentField)
elasticityDependentField.VariableLabelSet(iron.FieldVariableTypes.U,"ElasticityDependent")
elasticityDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"ElasticityTraction")
# Set the pressure to be nodally based and use the second mesh component
elasticityDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
elasticityDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
elasticityEquationsSet.DependentCreateFinish()

# Initialise elasticity dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
    elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
    elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
iron.Field.ComponentValuesInitialiseDP(
    elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,pInit)

# Create a field for the stress field. We will use this as the independent field for fitting to save a field copy. 
elasticityStressField = iron.Field()
elasticityStressField.CreateStart(elasticityStressFieldUserNumber,region)
elasticityStressField.TypeSet(iron.FieldTypes.GENERAL)
elasticityStressField.MeshDecompositionSet(decomposition)
elasticityStressField.GeometricFieldSet(geometricField)
elasticityStressField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
elasticityStressField.NumberOfVariablesSet(2)
elasticityStressField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.V])
elasticityStressField.VariableLabelSet(iron.FieldVariableTypes.U,"GaussStress")
elasticityStressField.VariableLabelSet(iron.FieldVariableTypes.V,"GaussWeight")
elasticityStressField.NumberOfComponentsSet(iron.FieldVariableTypes.U,6)
elasticityStressField.NumberOfComponentsSet(iron.FieldVariableTypes.V,6)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,5,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,6,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.V,1,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.V,2,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.V,3,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.V,4,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.V,5,1)
elasticityStressField.ComponentMeshComponentSet(iron.FieldVariableTypes.V,6,1)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.U,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.U,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.V,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.V,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.V,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.V,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.V,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.ComponentInterpolationSet(iron.FieldVariableTypes.V,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
elasticityStressField.CreateFinish()

# Create the derived equations set stress fields
elasticityEquationsSet.DerivedCreateStart(elasticityStressFieldUserNumber,elasticityStressField)
elasticityEquationsSet.DerivedVariableSet(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS,iron.FieldVariableTypes.U)
elasticityEquationsSet.DerivedCreateFinish()

# Create the material field
elasticityMaterialsField = iron.Field()
elasticityEquationsSet.MaterialsCreateStart(elasticityMaterialsFieldUserNumber,elasticityMaterialsField)
elasticityMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
elasticityMaterialsField.VariableLabelSet(iron.FieldVariableTypes.V,"Density")
elasticityEquationsSet.MaterialsCreateFinish()

# Set materials parameters
elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                     1,mooneyRivlin1)
elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                     2,mooneyRivlin2)
elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES, \
                                                     1,density)

# Create elasticity equations
elasticityEquations = iron.Equations()
elasticityEquationsSet.EquationsCreateStart(elasticityEquations)
elasticityEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
elasticityEquations.outputType = iron.EquationsOutputTypes.NONE
elasticityEquationsSet.EquationsCreateFinish()

# Create the fitting equations_set
fittingEquationsSetField = iron.Field()
fittingEquationsSet = iron.EquationsSet()
fittingEquationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                                    iron.EquationsSetTypes.GAUSS_FITTING_EQUATION,
                                    iron.EquationsSetSubtypes.GAUSS_POINT_FITTING,
                                    iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
fittingEquationsSet.CreateStart(fittingEquationsSetUserNumber,region,geometricField,
                                fittingEquationsSetSpecification,fittingEquationsSetFieldUserNumber,
                                fittingEquationsSetField)
fittingEquationsSet.CreateFinish()

# Create the fitting dependent field
fittingDependentField = iron.Field()
fittingEquationsSet.DependentCreateStart(fittingDependentFieldUserNumber,fittingDependentField)
fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.U,"FittedStress")
# Set the number of components to 2
fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,6)
fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,6)
# Set the field variables to be triquadratic Lagrange
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,5,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,6,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,5,1)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,6,1)
# Finish creating the fitting dependent field
fittingEquationsSet.DependentCreateFinish()

# Create the fitting independent field. Use the previously created elasticity stress field
fittingEquationsSet.IndependentCreateStart(elasticityStressFieldUserNumber,elasticityStressField)
# Finish creating the fitting independent field
fittingEquationsSet.IndependentCreateFinish()

# Initialise Gauss point weight field to 1.0
elasticityStressField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,1.0)
elasticityStressField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,2,1.0)
elasticityStressField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,3,1.0)
elasticityStressField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,4,1.0)
elasticityStressField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,5,1.0)
elasticityStressField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,6,1.0)

# Create material field (Sobolev parameters)
fittingMaterialField = iron.Field()
fittingEquationsSet.MaterialsCreateStart(fittingMaterialsFieldUserNumber,fittingMaterialField)
fittingMaterialField.VariableLabelSet(iron.FieldVariableTypes.U,"SmoothingParameters")
fittingEquationsSet.MaterialsCreateFinish()
# Set kappa and tau - Sobolev smoothing parameters
fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

# Create the fitting equations
fittingEquations = iron.Equations()
fittingEquationsSet.EquationsCreateStart(fittingEquations)
# Set the fitting equations sparsity type
fittingEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
# Set the fitting equations output type to none
fittingEquations.outputType = iron.EquationsOutputTypes.NONE
# Finish creating the fitting equations
fittingEquationsSet.EquationsCreateFinish()

# Define the elasticity problem
elasticityProblem = iron.Problem()
elasticityProblemSpecification = [iron.ProblemClasses.ELASTICITY,
                                  iron.ProblemTypes.FINITE_ELASTICITY,
                                  iron.ProblemSubtypes.STATIC_FINITE_ELASTICITY]
elasticityProblem.CreateStart(elasticityProblemUserNumber,iron.Context,elasticityProblemSpecification)
elasticityProblem.CreateFinish()

# Create the elasticity problem control loop
elasticityProblem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
elasticityProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
elasticityProblem.ControlLoopCreateFinish()

# Create elasticity problem solvers
elasticityNonLinearSolver = iron.Solver()
elasticityLinearSolver = iron.Solver()
elasticityProblem.SolversCreateStart()
elasticityProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,elasticityNonLinearSolver)
elasticityNonLinearSolver.outputType = iron.SolverOutputTypes.MONITOR
#elasticityNonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
#elasticityNonLinearSolver.outputType = iron.SolverOutputTypes.MATRIX
elasticityNonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
elasticityNonLinearSolver.NewtonAbsoluteToleranceSet(1e-14)
elasticityNonLinearSolver.NewtonSolutionToleranceSet(1e-14)
elasticityNonLinearSolver.NewtonRelativeToleranceSet(1e-14)
elasticityNonLinearSolver.NewtonLinearSolverGet(elasticityLinearSolver)
#elasticityLinearSolver.linearType = iron.LinearSolverTypes.DIRECT
elasticityProblem.SolversCreateFinish()

# Create elasticity solver equations and add elasticity equations set to solver equations
elasticitySolverEquations = iron.SolverEquations()
elasticityProblem.SolverEquationsCreateStart()
elasticityNonLinearSolver.SolverEquationsGet(elasticitySolverEquations)
elasticitySolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
elasticityEquationsSetIndex = elasticitySolverEquations.EquationsSetAdd(elasticityEquationsSet)
elasticityProblem.SolverEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
elasticityBoundaryConditions = iron.BoundaryConditions()
elasticitySolverEquations.BoundaryConditionsCreateStart(elasticityBoundaryConditions)

for widthNodeIdx in range(1,numberOfXNodes+1):
    for heightNodeIdx in range(1,numberOfYNodes+1):
        # Set left hand build in nodes ot no displacement
        nodeIdx=widthNodeIdx+(heightNodeIdx-1)*numberOfXNodes
        elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeIdx,1,
                                             iron.BoundaryConditionsTypes.FIXED,0.0)
        elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeIdx,2,
                                             iron.BoundaryConditionsTypes.FIXED,0.0)
        elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeIdx,3,
                                             iron.BoundaryConditionsTypes.FIXED,0.0)
        print("Wall node number = %d" % (nodeIdx))
    # Set downward force on right-hand edge
    nodeIdx=numberOfNodes-widthNodeIdx+1
    elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeIdx,2,
                                         iron.BoundaryConditionsTypes.NEUMANN_POINT,force)
    print("Force node number = %d" % (nodeIdx))
# Set reference pressure
elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,numberOfNodes,4,
                                     iron.BoundaryConditionsTypes.FIXED,pRef)

elasticitySolverEquations.BoundaryConditionsCreateFinish()

# Create fitting problem
fittingProblemSpecification = [iron.ProblemClasses.FITTING,
                               iron.ProblemTypes.DATA_FITTING,
                               iron.ProblemSubtypes.STATIC_FITTING]
fittingProblem = iron.Problem()
fittingProblem.CreateStart(fittingProblemUserNumber,iron.Context,fittingProblemSpecification)
fittingProblem.CreateFinish()

# Create control loops
fittingProblem.ControlLoopCreateStart()
fittingProblem.ControlLoopCreateFinish()

# Create problem solver
fittingSolver = iron.Solver()
fittingProblem.SolversCreateStart()
fittingProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,fittingSolver)
fittingSolver.outputType = iron.SolverOutputTypes.PROGRESS
fittingProblem.SolversCreateFinish()

# Create fitting solver equations and add fitting equations set to solver equations
fittingSolverEquations = iron.SolverEquations()
fittingProblem.SolverEquationsCreateStart()
# Get the solver equations
fittingSolver.SolverEquationsGet(fittingSolverEquations)
fittingSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
fittingEquationsSetIndex = fittingSolverEquations.EquationsSetAdd(fittingEquationsSet)
fittingProblem.SolverEquationsCreateFinish()

# Prescribe boundary conditions for the fitting problem
fittingBoundaryConditions = iron.BoundaryConditions()
fittingSolverEquations.BoundaryConditionsCreateStart(fittingBoundaryConditions)
fittingSolverEquations.BoundaryConditionsCreateFinish()

# Solve the elasticity problem
elasticityProblem.Solve()

# Calculate the stress field
elasticityEquationsSet.DerivedVariableCalculate(iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS)
                        
# Solve the fitting problem
fittingProblem.Solve()
                    
# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("Cantilever","FORTRAN")
fields.ElementsExport("Cantilever","FORTRAN")
fields.Finalise()

