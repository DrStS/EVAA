*** Settings ***
Library     CompareCsv
Library     OperatingSystem

*** Variables ***
${ABSTOLERANCE}     1e-12
${RELTOLERANCE}     1e-6

*** Test Cases ***
Abaqus
    ${ref}              Set Variable    ${CURDIR}/references/abaqus/Abaqus_11DofHomDiBC.dat
    ${evaa_input}       Set Variable    ${CURDIR}/../inputFiles/inputFileLinear.xml

    
    ${Res}              Run             ${CURDIR}/../build/bin/EVAA -i ${evaa_input} > ${CURDIR}/${TEST NAME}.dat

    Compare End State   ${ref}          ${CURDIR}/${TEST NAME}.dat    ${ABSTOLERANCE}     ${RELTOLERANCE}	

