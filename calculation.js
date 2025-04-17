export const standardGaussPoints = [
    { xi: -1/Math.sqrt(3), eta: -1/Math.sqrt(3) },
    { xi: 1/Math.sqrt(3), eta: -1/Math.sqrt(3) },
    { xi: 1/Math.sqrt(3), eta: 1/Math.sqrt(3) },
    { xi: -1/Math.sqrt(3), eta: 1/Math.sqrt(3) }
];

// Material properties
const E = 210000; // Young's modulus
const nu = 0.3;   // Poisson's ratio
const t = 1.0;    // Thickness


// Shape Function Calculations
export function calculateShapeFunctions(xi, eta) {
    return [
        0.25 * (1 - xi) * (1 - eta),  // N1
        0.25 * (1 + xi) * (1 - eta),  // N2
        0.25 * (1 + xi) * (1 + eta),  // N3
        0.25 * (1 - xi) * (1 + eta)   // N4
    ];
}

export function calculateDNdXi(xi, eta) {
    return [
        [-0.25 * (1 - eta), -0.25 * (1 - xi)],  // dN1/dξ, dN1/dη
        [0.25 * (1 - eta), -0.25 * (1 + xi)],   // dN2/dξ, dN2/dη
        [0.25 * (1 + eta), 0.25 * (1 + xi)],    // dN3/dξ, dN3/dη
        [-0.25 * (1 + eta), 0.25 * (1 - xi)]    // dN4/dξ, dN4/dη
    ];
}

export function calculateJacobiMatrix(dNdXi, nodes) {
    const J = [[0, 0], [0, 0]];

    for(let i = 0; i < 4; i++) {
        J[0][0] += dNdXi[i][0] * nodes[i].x;  // dx/dξ
        J[0][1] += dNdXi[i][1] * nodes[i].x;  // dx/dη
        J[1][0] += dNdXi[i][0] * nodes[i].y;  // dy/dξ
        J[1][1] += dNdXi[i][1] * nodes[i].y;  // dy/dη
    }

    return J;
}

export function calculateDeterminant(J) {
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

export function calculateInverseJacobi(J, detJ) {
    return [
        [J[1][1] / detJ, -J[0][1] / detJ],
        [-J[1][0] / detJ, J[0][0] / detJ]
    ];
}

export function calculateDNdX(dNdXi, Jinv) {
    const dNdX = [];
    for(let i = 0; i < 4; i++) {
        const dNdx = dNdXi[i][0] * Jinv[0][0] + dNdXi[i][1] * Jinv[0][1];
        const dNdy = dNdXi[i][0] * Jinv[1][0] + dNdXi[i][1] * Jinv[1][1];
        dNdX.push([dNdx, dNdy]);
    }
    return dNdX;
}

export function calculateBMatrix(dNdX) {
    const nN = dNdX.length; // number of nodes
    // Initialize B matrix with zeros (3 rows, 2*nN columns)
    const B = Array(3).fill().map(() => Array(2 * nN).fill(0));

    // Fill B matrix according to the MATLAB pattern:
    // Row 1: dN/dx terms in odd columns
    // Row 2: dN/dy terms in even columns
    // Row 3: dN/dy terms in odd columns and dN/dx terms in even columns

    for(let i = 0; i < nN; i++) {
        // First row: dN/dx terms
        B[0][2*i] = dNdX[i][0];     // dN/dx in odd columns

        // Second row: dN/dy terms
        B[1][2*i + 1] = dNdX[i][1]; // dN/dy in even columns

        // Third row: mixed terms
        B[2][2*i] = dNdX[i][1];     // dN/dy in odd columns
        B[2][2*i + 1] = dNdX[i][0]; // dN/dx in even columns
    }

    return B;
}

// Plane stress D matrix
export function calculateDMatrix() {
    const factor = E / (1 - nu * nu);
    return [
        [factor, factor * nu, 0],
        [factor * nu, factor, 0],
        [0, 0, factor * (1 - nu) / 2]
    ];
}

export function calculateStiffnessMatrix(nodes) {
    //todo is the weight from gauss in here?
    const D = calculateDMatrix();
    const nNodes = nodes.length;
    const matrixSize = 2 * nNodes; // 2 DOF per node

    // Initialize matrix with zeros
    const K = Array(matrixSize).fill().map(() => Array(matrixSize).fill(0));
    const gaussWeight = 1; // Weight for 2x2 Gauss quadrature

    const contributions = [];

    standardGaussPoints.forEach((gaussPoint) => {
        const {xi, eta} = gaussPoint;

        // Calculate matrices for this Gauss point
        const dNdXi = calculateDNdXi(xi, eta);
        const J = calculateJacobiMatrix(dNdXi, nodes);
        const detJ = calculateDeterminant(J);
        const Jinv = calculateInverseJacobi(J, detJ);
        const dNdX = calculateDNdX(dNdXi, Jinv);
        const B = calculateBMatrix(dNdX);

        // Calculate B transpose
        const BT = transposeMatrix(B);

        // Calculate B^T * D
        const BTD = multiplyMatrices(BT, D);

        // Calculate (B^T * D) * B
        const BTDB = multiplyMatrices(BTD, B);

        // Multiply by determinant and Gauss weight
        const contribution = BTDB.map(row =>
            row.map(value => value * detJ * gaussWeight)
        );

        // Add contribution to global stiffness matrix
        for(let i = 0; i < matrixSize; i++) {
            for(let j = 0; j < matrixSize; j++) {
                K[i][j] += contribution[i][j];
            }
        }

        contributions.push({
            gaussPoint: {xi, eta},
            contribution: contribution,
            detJ: detJ
        });
    });

    return {K, contributions};
}

export function reducedStiffnessMatrix(nodes) {
    const K = calculateStiffnessMatrix(nodes).K;

    // 2. Apply boundary conditions by modifying K
    const fixedDOFs = [
        6, 7,    // Node 4 (both x and y fixed)
        5        // Node 3 (y fixed)
    ];

    // 5. Reduce matrices by removing fixed DOFs
    const freeDOFs = [];
    for(let i = 0; i < 8; i++) {
        if(!fixedDOFs.includes(i)) freeDOFs.push(i);
    }

    // 6. Create reduced matrices
    const Kr = [];
    for(let i of freeDOFs) {
        const row = [];
        for(let j of freeDOFs) {
            row.push(K[i][j]);
        }
        Kr.push(row);
    }
    return Kr;
}

export function calculateDisplacements(nodes, forceArrow) {
    const fixedDOFs = [
        6, 7,    // Node 4 (both x and y fixed)
        5        // Node 3 (y fixed)
    ];

    const Fx = forceArrow.length * Math.cos(forceArrow.angle);
    const Fy = forceArrow.length * Math.sin(forceArrow.angle);

    // 4. Create force vector (8x1)
    const F = new Array(8).fill(0);
    F[2] = Fx;  // Node 2 x-force
    F[3] = Fy;  // Node 2 y-force

    // 5. Reduce matrices by removing fixed DOFs
    const freeDOFs = [];
    for(let i = 0; i < 8; i++) {
        if(!fixedDOFs.includes(i)) freeDOFs.push(i);
    }

    // 6. Create reduced matrices
    const Kr = reducedStiffnessMatrix(nodes);
    const Fr = [];
    for(let i of freeDOFs) {
        Fr.push(F[i]);
    }
    // 7. Solve Kr * Ur = Fr for Ur using Gaussian elimination
    const Ur = solveSystem(Kr, Fr);

    // 8. Create full displacement vector
    const U = new Array(8).fill(0);
    freeDOFs.forEach((dof, i) => {
        U[dof] = Ur[i];
    });

    // 9. Return displacements per node
    return [
        {x: U[0], y: U[1]},  // Node 1
        {x: U[2], y: U[3]},  // Node 2
        {x: U[4], y: U[5]},  // Node 3
        {x: U[6], y: U[7]}   // Node 4
    ];
}

// Helper to format matrices/vectors for logging (optional but nice)
function formatForLog(matrix) {
    // Limits the number of decimal places for cleaner output
    const formatNumber = (num) => (typeof num === 'number' ? num.toFixed(4) : num);

    if (!Array.isArray(matrix)) return matrix; // Not an array

    if (Array.isArray(matrix[0])) {
        // It's a 2D array (matrix)
        return matrix.map(row =>
            '[' + row.map(formatNumber).join(', ') + ']'
        ).join('\n  ');
    } else {
        // It's a 1D array (vector) - display as column for consistency if needed
        return matrix.map(el => `[${formatNumber(el)}]`).join('\n  ');
        // Or simply as a row: '[' + matrix.map(formatNumber).join(', ') + ']'
    }
}


export function calculateStrainAndStress(nodes, displacements) {
    const D = calculateDMatrix();
    const results = [];

    standardGaussPoints.forEach((gaussPoint, index) => {
        const { xi, eta } = gaussPoint;

        // Calculate matrices for this Gauss point
        const dNdXi = calculateDNdXi(xi, eta);
        const J = calculateJacobiMatrix(dNdXi, nodes);
        const detJ = calculateDeterminant(J);
        if (detJ <= 0) {
            console.error(`[GP ${index + 1}] Invalid Jacobian determinant: ${detJ}. Skipping point.`);
            // Handle error appropriately - maybe push null or throw error
            // For now, just log and continue to avoid crashing if possible
             results.push({
                gaussPoint: {xi, eta},
                error: `Invalid Jacobian determinant: ${detJ}`
            });
            return; // Skip this Gauss point
        }
        const Jinv = calculateInverseJacobi(J, detJ);
        const dNdX = calculateDNdX(dNdXi, Jinv);
        const B = calculateBMatrix(dNdX);

        // Create displacement vector U (as a column vector/matrix)
        const U_vector = [];
        displacements.forEach(d => {
            U_vector.push([d.x]); // x displacement
            U_vector.push([d.y]); // y displacement
        });

        // Calculate strain: ε = B * U
        const strain = multiplyMatrices(B, U_vector); // strain will be [[εxx], [εyy], [γxy]]

        // Calculate stress: σ = D * ε
        const stress = multiplyMatrices(D, strain);

        // Create result object with stress components
        const stressComponents = {
            σxx: stress[0][0],
            σyy: stress[1][0],
            τxy: stress[2][0]
        };

        // Calculate von Mises stress
        const vonMises = calculateVonMisesStress(stressComponents);

        results.push({
            gaussPoint: { xi, eta },
            strain: {
                εxx: strain[0][0],
                εyy: strain[1][0],
                γxy: strain[2][0] // Note: This is engineering shear strain (gamma_xy)
            },
            stress: stressComponents,
            vonMises: vonMises // Add von Mises stress to results
        });
    });

    return results;
}

// Calculate von Mises stress
export function calculateVonMisesStress(stress) {
    const {σxx, σyy, τxy} = stress;
    return Math.sqrt(σxx*σxx - σxx*σyy + σyy*σyy + 3*τxy*τxy);
}



// Helper functions for matrix operations
function transposeMatrix(matrix) {
    return matrix[0].map((_, colIndex) =>
        matrix.map(row => row[colIndex])
    );
}

function multiplyMatrices(a, b) {
    // Input validation
    if (!a || !b || !Array.isArray(a) || !Array.isArray(b)) {
        throw new Error('Invalid input: Both arguments must be arrays');
    }

    if (a.length === 0 || b.length === 0) {
        throw new Error('Invalid input: Empty matrices');
    }

    if (!Array.isArray(b[0])) {
        // If b is a vector, convert it to a column matrix
        b = b.map(value => [value]);
    }

    // Check if matrices can be multiplied
    const aColumns = a[0].length;
    const bRows = b.length;

    if (aColumns !== bRows) {
        throw new Error(`Matrix dimensions don't match: ${a.length}x${aColumns} and ${bRows}x${b[0].length}`);
    }

    // Perform multiplication
    return a.map(row => {
        return b[0].map((_, j) => {
            return row.reduce((sum, element, i) => {
                return sum + element * (b[i][j] || 0);
            }, 0);
        });
    });
}

function solveSystem(A, b) {
    const n = A.length;
    const x = new Array(n).fill(0);

    // Forward elimination
    for(let i = 0; i < n; i++) {
        for(let j = i + 1; j < n; j++) {
            const factor = A[j][i] / A[i][i];
            for(let k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // Back substitution
    for(let i = n - 1; i >= 0; i--) {
        let sum = 0;
        for(let j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

// Helper for matrix-vector multiplication (Ax = b)
function multiplyMatrixVector(matrix, vector) {
    const result = [0, 0, 0, 0];
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

// Extrapolation function using the precomputed inverse
export function extrapolateGaussToNodes(gaussResults, componentKey) {
    // componentKey is like 'Sxx', 'Eyy', etc.
    const gaussValues = [];

    if (!gaussResults || gaussResults.length !== 4) {
        console.error("Invalid gaussResults for extrapolation");
        return [0, 0, 0, 0]; // Return zero array on error
    }

    // Determine where the value lives (stress, strain, vonMises)
    let path;
    if (componentKey === 'Svm') {
         path = ['vonMises'];
    } else if (componentKey.startsWith('S')) {
        path = ['stress', componentKey === 'Sxx' ? 'σxx' : componentKey === 'Syy' ? 'σyy' : 'τxy'];
    } else if (componentKey.startsWith('E')) {
        path = ['strain', componentKey === 'Exx' ? 'εxx' : componentKey === 'Eyy' ? 'εyy' : 'γxy'];
    } else {
        console.error("Unknown component key for extrapolation:", componentKey);
        return [0, 0, 0, 0];
    }

    try {
        for (let i = 0; i < 4; i++) {
            let value = gaussResults[i];
            for (const p of path) {
                value = value[p];
            }
             if (typeof value !== 'number' || !isFinite(value)) {
                 console.warn(`Invalid value found during extrapolation for ${componentKey} at GP ${i}:`, value);
                 value = 0; // Assign 0 if invalid
             }
            gaussValues.push(value);
        }
    } catch (e) {
        console.error(`Error accessing component ${componentKey} in gaussResults:`, e);
        return [0, 0, 0, 0];
    }


    // Perform extrapolation: NodalValues = N_inv * GaussValues
    const N_inv = [
      [ 1.8660254, -0.5      ,  0.1339746, -0.5      ],
      [-0.5      ,  1.8660254, -0.5      ,  0.1339746],
      [ 0.1339746, -0.5      ,  1.8660254, -0.5      ],
      [-0.5      ,  0.1339746, -0.5      ,  1.8660254]
    ];
    const nodalValues = multiplyMatrixVector(N_inv, gaussValues);
    return nodalValues;
}