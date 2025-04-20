// --- Gauss Quadrature Definitions ---
export const gaussQuadratureRules = {
    '1x1': [
        { xi: 0, eta: 0, weight: 4 } // For [-1, 1] x [-1, 1], the weight is 2 * 2 = 4
    ],
    '2x2': [
        { xi: -1/Math.sqrt(3), eta: -1/Math.sqrt(3), weight: 1 },
        { xi:  1/Math.sqrt(3), eta: -1/Math.sqrt(3), weight: 1 },
        { xi:  1/Math.sqrt(3), eta:  1/Math.sqrt(3), weight: 1 },
        { xi: -1/Math.sqrt(3), eta:  1/Math.sqrt(3), weight: 1 }
    ]
};

// --- Material Properties (Example - should be passed in or configured) ---
const defaultMaterialProps = {
    E: 210000, // Young's modulus
    nu: 0.3,   // Poisson's ratio
    t: 1.0    // Thickness (assuming plane stress/strain context needs it)
};

// --- Shape Function Calculations ---
export function calculateShapeFunctions(xi, eta) {
    return [
        0.25 * (1 - xi) * (1 - eta),  // N1
        0.25 * (1 + xi) * (1 - eta),  // N2
        0.25 * (1 + xi) * (1 + eta),  // N3
        0.25 * (1 - xi) * (1 + eta)   // N4
    ];
}

// Calculate Shape Function Derivatives w.r.t. local coords (ksi, eta)
// Returns: [[dN1/dxi, dN1/deta], [dN2/dxi, dN2/deta], ...] (4x2)
export function calculateDNdXi(xi, eta) {
    return [
        [-0.25 * (1 - eta), -0.25 * (1 - xi)],  // dN1/dξ, dN1/dη
        [ 0.25 * (1 - eta), -0.25 * (1 + xi)],  // dN2/dξ, dN2/dη
        [ 0.25 * (1 + eta),  0.25 * (1 + xi)],  // dN3/dξ, dN3/dη
        [-0.25 * (1 + eta),  0.25 * (1 - xi)]   // dN4/dξ, dN4/dη
    ];
}

// Calculate STANDARD Jacobian Matrix J = [[dx/dxi, dx/deta], [dy/dxi, dy/deta]] (2x2)
// nodes: Array of node coordinate objects [{x, y}, ...]
// dNdXi: Array of shape function derivatives [[dN1/dxi, dN1/deta], ...] (4x2)
export function calculateJacobiMatrix(dNdXi, nodes) {
    const J = [[0, 0], [0, 0]];
    if (nodes.length !== dNdXi.length) {
        throw new Error("Number of nodes must match the number of shape functions derivatives.");
    }
    const nNodes = nodes.length;
    for(let i = 0; i < nNodes; i++) {
        J[0][0] += dNdXi[i][0] * nodes[i].x;  // dx/dξ = sum(dNi/dξ * xi)
        J[0][1] += dNdXi[i][1] * nodes[i].x;  // dx/dη = sum(dNi/dη * xi)
        J[1][0] += dNdXi[i][0] * nodes[i].y;  // dy/dξ = sum(dNi/dξ * yi)
        J[1][1] += dNdXi[i][1] * nodes[i].y;  // dy/dη = sum(dNi/dη * yi)
    }
    return J;
}

// Calculate Determinant of a 2x2 Matrix
export function calculateDeterminant(J) {
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

// Calculate Inverse of a 2x2 Jacobian Matrix
export function calculateInverseJacobi(J, detJ) {
    if (Math.abs(detJ) < 1e-12) { // Check for near-zero determinant
        console.warn("Near-zero Jacobian determinant encountered:", detJ);
        // Return an identity matrix or throw an error, depending on desired handling
        // For now, returning a matrix that might lead to NaNs downstream to flag the issue
        return [ [Infinity, Infinity], [Infinity, Infinity] ];
    }
    const invDetJ = 1.0 / detJ;
    return [
        [ J[1][1] * invDetJ, -J[0][1] * invDetJ], // dxi/dx,  dxi/dy
        [-J[1][0] * invDetJ,  J[0][0] * invDetJ]  // deta/dx, deta/dy
    ];
}

// Calculate Shape Function Derivatives w.r.t. global coords (x, y)
// Returns: [[dN1/dx, dN1/dy], [dN2/dx, dN2/dy], ...] (4x2)
export function calculateDNdX(dNdXi, Jinv) {
    const dNdX = [];
    const nNodes = dNdXi.length;
    for(let i = 0; i < nNodes; i++) {
        const dNdxi = dNdXi[i][0]; // dNi/dξ
        const dNdeta = dNdXi[i][1]; // dNi/dη
        // dNi/dx = dNi/dξ * dξ/dx + dNi/dη * dη/dx
        const dNdx = dNdxi * Jinv[0][0] + dNdeta * Jinv[1][0];
        // dNi/dy = dNi/dξ * dξ/dy + dNi/dη * dη/dy
        const dNdy = dNdxi * Jinv[0][1] + dNdeta * Jinv[1][1];
        dNdX.push([dNdx, dNdy]);
    }
    return dNdX;
}

// Calculate Strain-Displacement Matrix B (3x8 for 4-node element)
export function calculateBMatrix(dNdX) {
    const nNodes = dNdX.length; // Should be 4 for Q4
    const nDofNode = 2;
    const nDofTotal = nNodes * nDofNode; // 8
    const B = Array(3).fill().map(() => Array(nDofTotal).fill(0));

    for(let i = 0; i < nNodes; i++) {
        const dNdx = dNdX[i][0];
        const dNdy = dNdX[i][1];
        const col1 = 2 * i;     // Column index for u_i displacement effect
        const col2 = 2 * i + 1; // Column index for v_i displacement effect

        // Row 1: ε_xx = sum(dNi/dx * ui)
        B[0][col1] = dNdx;

        // Row 2: ε_yy = sum(dNi/dy * vi)
        B[1][col2] = dNdy;

        // Row 3: γ_xy = sum(dNi/dy * ui + dNi/dx * vi)
        B[2][col1] = dNdy;
        B[2][col2] = dNdx;
    }
    return B;
}

// Calculate Plane Stress Constitutive Matrix D
export function calculateDMatrix(materialProps = defaultMaterialProps) {
    const { E, nu } = materialProps;
    if (typeof E !== 'number' || typeof nu !== 'number' || !isFinite(E) || !isFinite(nu)) {
        console.error("Invalid E or nu in material properties:", materialProps);
        throw new Error("Invalid material properties for D matrix calculation.");
    }
    // Calculate G (Shear Modulus) here as well, used in D[2][2]
    const G = E / (2.0 * (1.0 + nu));
    const denom = 1.0 - nu * nu;
    if (Math.abs(denom) < 1e-14) {
         throw new Error("Cannot calculate D matrix: 1 - nu^2 is close to zero.");
    }
    const factor = E / denom;

    return [
        [factor,        factor * nu,   0 ],
        [factor * nu,   factor,        0 ],
        [0,             0,             G ] // Use G directly = E / (2*(1+nu))
    ];
}

/**
 * Calculates the stiffness matrix for a 2D quadrilateral element.
 * @param {Array<object>} nodes - Array of node coordinates [{x, y}, {x, y}, {x, y}, {x, y}]
 * @param {object} materialProps - Object containing material properties (E, nu, t). Uses defaults if not provided.
 * @param {string} quadratureRuleName - The name of the quadrature rule ('1x1' or '2x2'). Defaults to '2x2'.
 * @param {number} hourglassAlpha - Stabilization parameter (e.g., 0.01). Only used for '1x1'. Defaults to 0.01.
 * @returns {object} Object containing the stiffness matrix K, contributions per Gauss point, and hourglass details (if applicable).
 */
export function calculateStiffnessMatrix(nodes, materialProps = defaultMaterialProps, quadratureRuleName = '2x2', hourglassAlpha = 0.05) { // Default alpha changed slightly
    // --- Input Validation and Setup ---
    if (!nodes || nodes.length !== 4) {
        throw new Error("Stiffness matrix calculation currently supports only 4-node quadrilateral elements.");
    }
    const selectedGaussPoints = gaussQuadratureRules[quadratureRuleName];
    if (!selectedGaussPoints) {
        throw new Error(`Invalid quadrature rule name: ${quadratureRuleName}. Use '1x1' or '2x2'.`);
    }
    if (!materialProps || typeof materialProps.E !== 'number' || typeof materialProps.nu !== 'number') {
         console.warn("Invalid or missing material properties. Using default values.");
         materialProps = { ...defaultMaterialProps }; // Use a copy of defaults
    } else {
        // Ensure E and nu are valid numbers before proceeding
        if (!isFinite(materialProps.E) || !isFinite(materialProps.nu)) {
             console.error("Non-finite E or nu provided:", materialProps);
             throw new Error("Non-finite material properties (E, nu) provided.");
        }
        // Ensure materialProps is a fresh object if defaults were merged or passed directly
        materialProps = { ...materialProps };
    }

    // Use thickness 't' from materialProps, default to 1.0 if not provided
    const thickness = typeof materialProps.t === 'number' && isFinite(materialProps.t) ? materialProps.t : 1.0;
    materialProps.t = thickness; // Ensure thickness is stored back if defaulted

    console.log(`Calculating Stiffness Matrix using ${quadratureRuleName} Gauss Quadrature (${selectedGaussPoints.length} points)`);
    if (quadratureRuleName === '1x1' && hourglassAlpha > 0) {
        // Check if alpha is reasonable
        if (!isFinite(hourglassAlpha) || hourglassAlpha < 0) {
             console.warn(`Invalid hourglassAlpha (${hourglassAlpha}). Setting to 0 (disabled).`);
             hourglassAlpha = 0;
        } else {
             console.log(` -> Applying Hourglass Stabilization with alpha = ${hourglassAlpha}`);
        }
    }

    const D = calculateDMatrix(materialProps); // Constitutive matrix (Plane Stress)

    const nNodes = 4;
    const nDofNode = 2;
    const matrixSize = nNodes * nDofNode; // 8

    // Initialize element stiffness matrix
    const K = Array(matrixSize).fill(0).map(() => Array(matrixSize).fill(0));
    const contributions = []; // Stores details for each Gauss point contribution

    // --- Standard Stiffness Integration Loop ---
    selectedGaussPoints.forEach((gaussPoint, index) => {
        const { xi, eta, weight } = gaussPoint;

        try {
            // 1. Derivatives w.r.t. local coordinates (ksi, eta)
            const dNdXi = calculateDNdXi(xi, eta); // 4x2

            // 2. Jacobian Matrix (Standard Definition)
            const J = calculateJacobiMatrix(dNdXi, nodes); // 2x2
            const detJ = calculateDeterminant(J);

            if (detJ <= 1e-12) { // Check for zero or negative determinant
                console.error(`[GP ${index+1}] Warning: Non-positive or near-zero Jacobian determinant (${detJ.toFixed(4)}) at (xi=${xi.toFixed(4)}, eta=${eta.toFixed(4)}). Element might be distorted or invalid. Skipping contribution.`);
                contributions.push({
                    gaussPoint: { xi, eta, weight },
                    contribution: Array(matrixSize).fill(0).map(() => Array(matrixSize).fill(0)),
                    detJ: detJ,
                    error: `Non-positive/near-zero Jacobian determinant: ${detJ}`
                });
                return; // Skip calculations for this Gauss point
            }

            // 3. Inverse Jacobian
            const Jinv = calculateInverseJacobi(J, detJ); // 2x2
             // Check if inverse calculation failed (returned Infinities)
             if (!Jinv.every(row => row.every(isFinite))) {
                  console.error(`[GP ${index+1}] Failed to compute inverse Jacobian (Jacobian likely singular) at (xi=${xi.toFixed(4)}, eta=${eta.toFixed(4)}). Skipping contribution.`);
                   contributions.push({
                        gaussPoint: { xi, eta, weight },
                        contribution: Array(matrixSize).fill(0).map(() => Array(matrixSize).fill(0)),
                        detJ: detJ,
                        error: `Failed to compute inverse Jacobian`
                    });
                    return;
             }

            // 4. Derivatives w.r.t. global coordinates (x, y)
            const dNdX = calculateDNdX(dNdXi, Jinv); // 4x2

            // 5. B Matrix (Strain-Displacement)
            const B = calculateBMatrix(dNdX); // 3x8

            // --- Calculate Standard Stiffness Contribution: K_gp = B^T * D * B * detJ * weight * thickness ---
            const BT = transposeMatrix(B);            // 8x3
            const BTD = multiplyMatrices(BT, D);      // 8x3
            const BTDB = multiplyMatrices(BTD, B);    // 8x8

            const integrationFactor = detJ * weight * thickness; // Combined factor dV

            const standardContribution = BTDB.map(row =>
                row.map(value => value * integrationFactor)
            );

            // --- Accumulate Contribution ---
            let nanEncountered = false;
            for (let i = 0; i < matrixSize; i++) {
                for (let j = 0; j < matrixSize; j++) {
                    if (isNaN(standardContribution[i][j]) || !isFinite(standardContribution[i][j])) {
                        console.error(`[GP ${index+1}] NaN or Infinity encountered in standard stiffness contribution at K[${i}][${j}] (xi=${xi.toFixed(4)}, eta=${eta.toFixed(4)}). Setting to 0.`);
                        standardContribution[i][j] = 0;
                        nanEncountered = true;
                    }
                    K[i][j] += standardContribution[i][j];
                }
            }

            // --- Store details for this Gauss point ---
            contributions.push({
                gaussPoint: { xi, eta, weight },
                contribution: standardContribution,
                detJ: detJ,
                B: B, // Optional: Store for debugging
                dNdX: dNdX, // Optional: Store for debugging
                error: nanEncountered ? "NaN/Infinity in calculation" : null
            });

        } catch (error) {
            console.error(`[GP ${index+1}] Error during calculation at (xi=${xi.toFixed(4)}, eta=${eta.toFixed(4)}):`, error);
            contributions.push({
                gaussPoint: { xi, eta, weight },
                contribution: Array(matrixSize).fill(0).map(() => Array(matrixSize).fill(0)),
                error: `Calculation error: ${error.message}`
            });
        }
    }); // End of loop over Gauss points


    // --- Hourglass Stabilization (Only for 1x1 rule and alpha > 0) ---
    let K_hg = null; // Initialize Hourglass stiffness contribution matrix
    let hgDetails = null; // Initialize details object

    if (quadratureRuleName === '1x1' && hourglassAlpha > 0 && contributions.length > 0 && !contributions[0].error) {
        // Proceed only if 1x1 rule, alpha > 0, and the standard 1x1 calc didn't fail

        console.log("   Applying Hourglass Stabilization Calculation...");
        K_hg = Array(matrixSize).fill(0).map(() => Array(matrixSize).fill(0)); // 8x8 matrix for HG stiffness

        try {
            // --- Calculations are done at the element center (xi=0, eta=0) ---
            const xi_center = 0.0;
            const eta_center = 0.0;
            const { E, nu } = materialProps; // Get E, nu from potentially updated props

            // Check E, nu again specifically for HG calc
            if (typeof E !== 'number' || typeof nu !== 'number' || !isFinite(E) || !isFinite(nu)) {
                 throw new Error("Invalid material properties (E, nu) for hourglass calculation.");
            }


            // 1. Calculate derivatives and Jacobian AT THE CENTER
            const dNdXi_center = calculateDNdXi(xi_center, eta_center); // 4x2
            const J_center = calculateJacobiMatrix(dNdXi_center, nodes); // 2x2
            const detJ_center = calculateDeterminant(J_center);

            if (detJ_center <= 1e-12) {
                console.error(`   Hourglass Control Warning: Non-positive/near-zero Jacobian determinant (${detJ_center.toFixed(4)}) at element center. Skipping stabilization.`);
                 hgDetails = { error: `Non-positive/near-zero Jacobian at center: ${detJ_center}` };
            } else {
                const Jinv_center = calculateInverseJacobi(J_center, detJ_center); // 2x2
                // Check inverse calculation
                 if (!Jinv_center.every(row => row.every(isFinite))) {
                      throw new Error("Failed to compute inverse Jacobian at center for hourglass calculation.");
                 }
                const dNdX_center = calculateDNdX(dNdXi_center, Jinv_center); // 4x2 [[dN1/dx, dN1/dy], ...]

                // 2. Define the base hourglass vector (gamma vector pattern)
                const gamma_base = [1.0, -1.0, 1.0, -1.0]; // Corresponds to nodes 1, 2, 3, 4

                // 3. Calculate the geometric projection vector 'b' (often called 'h' or 'xi' elsewhere)
                //    b_x = sum(gamma_base[a] * x_a), b_y = sum(gamma_base[a] * y_a)
                const b_vector = [0.0, 0.0]; // [b_x, b_y]
                for (let a = 0; a < nNodes; a++) { // Loop over nodes
                    b_vector[0] += gamma_base[a] * nodes[a].x; // Sum gamma[a] * x[a]
                    b_vector[1] += gamma_base[a] * nodes[a].y; // Sum gamma[a] * y[a]
                }

                // 4. Orthogonalize the base gamma vector -> gamma_corrected (or gamma_hat)
                //    gamma_corrected[a] = gamma_base[a] - sum_{i=x,y} ( b_i * dNa/di )
                let gamma_corrected = [...gamma_base]; // Start with a copy

                for (let a = 0; a < nNodes; a++) { // Loop over nodes
                    let correction_term = 0.0;
                    const dNdx_a = dNdX_center[a][0]; // dNa/dx
                    const dNdy_a = dNdX_center[a][1]; // dNa/dy

                    // Correction term calculation: b dot grad(Na)
                    correction_term += b_vector[0] * dNdx_a; // b_x * dNa/dx
                    correction_term += b_vector[1] * dNdy_a; // b_y * dNa/dy

                    gamma_corrected[a] -= correction_term; // Apply correction
                }

                // 5. Calculate the hourglass stiffness contribution MAGNITUDE
                //    *** This is the corrected part ***
                //    Magnitude proportional to alpha * G * Area * thickness
                const G = E / (2.0 * (1.0 + nu)); // Shear Modulus
                const Area = detJ_center * 4.0; // Area = detJ * weight (weight=4 for 1x1 on [-1,1]x[-1,1])
                // Check for valid Area and G
                if (!isFinite(G) || !isFinite(Area)) {
                     throw new Error(`Invalid G (${G}) or Area (${Area}) calculated for hourglass stiffness.`);
                }
                const hgStiffnessMagnitude = hourglassAlpha * G * Area * thickness;
                // Removed: 'scale_factor' calculation
                // Removed: stabFactor = hourglassAlpha * scale_factor * Area * thickness;

                // 6. Assemble K_hg (8x8) using corrected Kronecker-product logic
                //    K_hg_submatrix = hgStiffnessMagnitude * gamma_corrected[a] * gamma_corrected[b] * Identity(2x2)
                for (let a = 0; a < nNodes; a++) { // Node index 'a' (row node)
                    for (let b = 0; b < nNodes; b++) { // Node index 'b' (column node)
                        const gamma_ab_term = gamma_corrected[a] * gamma_corrected[b];
                        const blockStiffness = hgStiffnessMagnitude * gamma_ab_term;

                        if (!isFinite(blockStiffness)) {
                            console.warn(`   Hourglass Warning: Non-finite block stiffness calculated for nodes (${a+1}, ${b+1}). Skipping block.`);
                            continue; // Skip this block if magnitude/gamma calculation failed
                        }

                        // Calculate global DOF indices for the top-left corner of the 2x2 block
                        let baseRow = nDofNode * a; // 2*a
                        let baseCol = nDofNode * b; // 2*b

                        // Add the 2x2 identity block scaled by blockStiffness
                        // Only add to diagonal terms of the 2x2 block (Identity matrix effect)
                        K_hg[baseRow + 0][baseCol + 0] += blockStiffness; // K_xx += stiffness * 1
                        K_hg[baseRow + 1][baseCol + 1] += blockStiffness; // K_yy += stiffness * 1
                        // Off-diagonal terms of the 2x2 block (K_xy, K_yx) are zero because Identity is diagonal
                        // K_hg[baseRow + 0][baseCol + 1] += blockStiffness * 0;
                        // K_hg[baseRow + 1][baseCol + 0] += blockStiffness * 0;
                    }
                }

                // 7. Add hourglass stiffness K_hg to the total stiffness K
                let nanInHg = false;
                for (let i = 0; i < matrixSize; i++) {
                    for (let j = 0; j < matrixSize; j++) {
                         // Final check before adding
                         if (isNaN(K_hg[i][j]) || !isFinite(K_hg[i][j])) {
                             console.error(`   Hourglass Control Error: NaN or Infinity value encountered in final K_hg[${i}][${j}]. Setting K_hg contribution to 0.`);
                             K_hg[i][j] = 0; // Prevent corrupting K
                             nanInHg = true;
                         }
                         K[i][j] += K_hg[i][j];
                    }
                }

                 // Store HG details for debugging/inspection
                 hgDetails = {
                    gamma_base: gamma_base,
                    b_vector: b_vector, // Renamed from h_vector
                    dNdX_center: dNdX_center,
                    gamma_corrected: gamma_corrected,
                    // scale_factor: undefined, // Removed
                    G: G,
                    Area: Area,
                    hgStiffnessMagnitude: hgStiffnessMagnitude, // Renamed from stabFactor
                    K_hg: K_hg, // The actual added HG stiffness matrix
                    error: nanInHg ? "NaN/Infinity in K_hg calculation" : null
                 };
                 console.log("   Hourglass stabilization applied successfully.");

            } // End if detJ_center > 0
        } catch (error) {
             console.error("   Error during hourglass stabilization calculation:", error);
             // Ensure K_hg remains null or zeroed if error occurs during calculation
             K_hg = Array(matrixSize).fill(0).map(() => Array(matrixSize).fill(0));
             hgDetails = { error: `Calculation error: ${error.message}` };
        }

    } // End of hourglass stabilization block


    // --- Return Results ---
    return { K, contributions, hgDetails };
}


// ========================================================================
// Other Functions (Assumed to be correct and used by the above)
// ========================================================================

// Placeholder/Example Helper functions for matrix operations
function transposeMatrix(matrix) {
    if (!matrix || matrix.length === 0 || !matrix[0]) return [];
    return matrix[0].map((_, colIndex) =>
        matrix.map(row => row[colIndex])
    );
}

function multiplyMatrices(A, B) {
    const rowsA = A.length;
    if (rowsA === 0) return [];
    const colsA = A[0].length;
    const rowsB = B.length;
    if (rowsB === 0) return new Array(rowsA).fill(null).map(() => []); // Handle B being empty
    const colsB = B[0].length;

    if (colsA !== rowsB) {
        throw new Error(`Cannot multiply matrices: dimensions mismatch ${rowsA}x${colsA} and ${rowsB}x${colsB}`);
    }

    const C = new Array(rowsA);
    for (let i = 0; i < rowsA; i++) {
        C[i] = new Array(colsB).fill(0);
        for (let j = 0; j < colsB; j++) {
            for (let k = 0; k < colsA; k++) {
                // Ensure values are numbers before multiplying
                const valA = typeof A[i][k] === 'number' ? A[i][k] : 0;
                const valB = typeof B[k][j] === 'number' ? B[k][j] : 0;
                C[i][j] += valA * valB;
            }
        }
    }
    return C;
}


// --- Functions below this line are assumed to exist and work correctly ---
// --- They are included for context but were not part of the refactoring request ---

export function reducedStiffnessMatrix(nodes, materialProps, quadratureRuleName = '2x2', hourglassAlpha = 0.01) {
    const { K } = calculateStiffnessMatrix(nodes, materialProps, quadratureRuleName, hourglassAlpha);

    // Example Fixed DOFs for a cantilever beam fixed at Node 4 (indices 6,7) and roller at Node 3 (index 5)
    const fixedDOFs = [
        5, // Node 3, y-DOF
        6, // Node 4, x-DOF
        7  // Node 4, y-DOF
    ];
    const matrixSize = K.length; // Should be 8

    const freeDOFs = [];
    for(let i = 0; i < matrixSize; i++) {
        if(!fixedDOFs.includes(i)) freeDOFs.push(i);
    }

    const Kr = [];
    for(let i = 0; i < freeDOFs.length; i++) {
        const row = [];
        for(let j = 0; j < freeDOFs.length; j++) {
            row.push(K[freeDOFs[i]][freeDOFs[j]]);
        }
        Kr.push(row);
    }
    return Kr;
}

// Placeholder solver - Use a robust library (like numeric.js or math.js) in practice!
function solveSystem(A, b) {
    console.warn("Using basic Gaussian elimination solver. Consider a robust library for stability.");
    const n = A.length;
    const Aug = A.map((row, i) => [...row, b[i]]); // Augmented matrix

    // Forward elimination
    for (let i = 0; i < n; i++) {
        // Find pivot
        let maxRow = i;
        for (let k = i + 1; k < n; k++) {
            if (Math.abs(Aug[k][i]) > Math.abs(Aug[maxRow][i])) {
                maxRow = k;
            }
        }
        // Swap rows
        [Aug[i], Aug[maxRow]] = [Aug[maxRow], Aug[i]];

        // Check for singularity
        if (Math.abs(Aug[i][i]) < 1e-12) {
             console.error("Matrix is singular or near-singular during forward elimination.");
            return new Array(n).fill(NaN); // Indicate failure
        }


        // Eliminate below
        for (let k = i + 1; k < n; k++) {
            const factor = Aug[k][i] / Aug[i][i];
            Aug[k][i] = 0; // Set element to 0 directly
            for (let j = i + 1; j <= n; j++) { // Include augmented part
                Aug[k][j] -= factor * Aug[i][j];
            }
        }
    }

    // Back substitution
    const x = new Array(n);
    for (let i = n - 1; i >= 0; i--) {
         if (Math.abs(Aug[i][i]) < 1e-12) {
             console.error("Matrix is singular or near-singular during back substitution.");
            return new Array(n).fill(NaN); // Indicate failure
        }
        let sum = 0;
        for (let j = i + 1; j < n; j++) {
            sum += Aug[i][j] * x[j];
        }
        x[i] = (Aug[i][n] - sum) / Aug[i][i];
    }

     // Final check for NaNs just in case
    if (x.some(isNaN)) {
        console.error("Solver resulted in NaN values.");
    }

    return x;
}


export function calculateDisplacements(nodes, forceArrow, materialProps, quadratureRuleName = '2x2', hourglassAlpha = 0.05) { // Match default alpha
    console.log("Calculating displacements using quadrature:", quadratureRuleName);
    // Example Fixed DOFs
     const fixedDOFs = [ 5, 6, 7 ]; // Fix N3y, N4x, N4y
     const matrixSize = 8;
     const nDofNode = 2;
     const nNodes = nodes.length;

     // --- Validate Inputs ---
     if (!nodes || nodes.length !== 4) {
         console.error("Invalid nodes input for displacement calculation.");
         return Array(4).fill({ x: 0, y: 0 }); // Return default shape
     }
      if (!materialProps || typeof materialProps.E !== 'number' || typeof materialProps.nu !== 'number' || !isFinite(materialProps.E) || !isFinite(materialProps.nu)) {
         console.warn("Invalid or missing material properties for displacement calc. Using defaults.");
         materialProps = { ...defaultMaterialProps };
     } else {
         materialProps = { ...materialProps }; // Ensure it's a copy
     }
      // Ensure thickness is present
      if (typeof materialProps.t !== 'number' || !isFinite(materialProps.t)) {
          materialProps.t = 1.0;
      }

    // --- Create Full Force Vector ---
    const F_full = new Array(matrixSize).fill(0);
    let appliedForceInfo = "Zero force applied.";
    if (forceArrow && typeof forceArrow.length === 'number' && typeof forceArrow.angle === 'number' && isFinite(forceArrow.length) && isFinite(forceArrow.angle)) {
        const Fx = forceArrow.length * Math.cos(forceArrow.angle);
        const Fy = forceArrow.length * Math.sin(forceArrow.angle);
         // Apply force to Node 2 (DOFs 2, 3) - Hardcoded for this example
         const forceNodeIndex = 1; // Node 2 (0-based)
         const forceDofX = nDofNode * forceNodeIndex;     // = 2
         const forceDofY = nDofNode * forceNodeIndex + 1; // = 3
         if (forceDofX < matrixSize && forceDofY < matrixSize) {
              F_full[forceDofX] = Fx;
              F_full[forceDofY] = Fy;
              appliedForceInfo = `Applied force (Fx=${Fx.toFixed(4)}, Fy=${Fy.toFixed(4)}) to Node ${forceNodeIndex + 1} (DOFs ${forceDofX}, ${forceDofY})`;
         } else {
              appliedForceInfo = "Warning: Force node index out of bounds. Zero force applied.";
              console.error(appliedForceInfo);
         }

    } else if (forceArrow) {
        appliedForceInfo = "Warning: Invalid forceArrow object (length/angle missing or non-finite). Zero force applied.";
        console.warn(appliedForceInfo, forceArrow);
    }
    console.log(appliedForceInfo);
    console.log("Full Force Vector (F_full):", F_full.map(f=>f.toFixed(4)));


    // --- Get Full Stiffness Matrix ---
    let K_full_obj;
    try {
         K_full_obj = calculateStiffnessMatrix(nodes, materialProps, quadratureRuleName, hourglassAlpha);
    } catch (error) {
        console.error("Error calculating stiffness matrix:", error);
         return Array(nNodes).fill({ x: 0, y: 0 }); // Return default shape on failure
    }
     const K_full = K_full_obj.K; // Extract the matrix

    // --- Get Reduced System ---
    const freeDOFs = [];
    const allDOFs = Array.from({length: matrixSize}, (_, i) => i);
    allDOFs.forEach(dof => {
        if (!fixedDOFs.includes(dof)) {
            freeDOFs.push(dof);
        }
    });

     if (freeDOFs.length === 0) {
        console.log("All DOFs are fixed. No system to solve. Returning zero displacements.");
        return Array(nNodes).fill({ x: 0, y: 0 });
    }

    // Extract reduced stiffness matrix Kr
    const Kr = freeDOFs.map(i =>
        freeDOFs.map(j => K_full[i][j])
    );
    // Extract reduced force vector Fr
    const Fr = freeDOFs.map(dof => F_full[dof]);

    console.log("Reduced Stiffness Matrix Kr:", Kr.map(row => row.map(v => v.toExponential(3))));
    console.log("Free DOFs:", freeDOFs);
    console.log("Reduced Force Vector Fr:", Fr.map(f => f.toFixed(4)));

    // --- Solve for Free Displacements ---
    console.log("Solving system Kr * Ur = Fr");
    let Ur;
    try {
        Ur = solveSystem(Kr, Fr); // Solve the reduced system Kr * Ur = Fr
    } catch (error) {
         console.error("Error solving the linear system Kr * Ur = Fr:", error);
         console.error("Kr:", Kr);
         console.error("Fr:", Fr);
         console.log("This often means Kr is singular due to insufficient boundary conditions or zero-energy modes.");
         return Array(nNodes).fill({ x: 0, y: 0 }); // Return zeros on solver failure
    }


    // --- Check Solver Result ---
    if (!Ur || Ur.some(u => isNaN(u) || !isFinite(u))) {
        console.error("Solver failed or produced non-finite displacements (NaN/Infinity).");
        console.log("Possible causes: Singular Kr (check fixedDOFs), numerical instability (check element geometry/material/HG), solver limitations.");
        // Return zero displacements to prevent downstream errors
        return Array(nNodes).fill({ x: 0, y: 0 });
    }
     console.log("Solved Free Displacements Ur:", Ur.map(u=>u.toExponential(4)));


    // --- Assemble Full Displacement Vector ---
    const U_full = new Array(matrixSize).fill(0); // Initialize full displacement vector with zeros
    freeDOFs.forEach((dof, i) => {
        if (i < Ur.length) { // Safety check
            U_full[dof] = Ur[i]; // Place calculated free displacements
        }
    });
    console.log("Full Displacement Vector U_full:", U_full.map(u=>u.toExponential(4)));

    // --- Format Output ---
    const displacements = [];
    for (let i = 0; i < nNodes; i++) {
        const dx = U_full[nDofNode * i];
        const dy = U_full[nDofNode * i + 1];
        displacements.push({ x: dx, y: dy });
    }

    return displacements;
}

export function calculateStrainAndStress(nodes, displacements, materialProps, quadratureRuleName = '2x2') {
    // ... (calculateStrainAndStress implementation - seems okay from snippet) ...
    // Ensure it correctly uses materialProps for D matrix
     const selectedGaussPoints = gaussQuadratureRules[quadratureRuleName];

    if (!selectedGaussPoints) {
        throw new Error(`Invalid quadrature rule name: ${quadratureRuleName}. Use '1x1' or '2x2'.`);
    }
     if (!materialProps) materialProps = defaultMaterialProps; // Use defaults if needed

     console.log(`Calculating Strain/Stress using ${quadratureRuleName} Gauss Quadrature (${selectedGaussPoints.length} points)`);

    const D = calculateDMatrix(materialProps); // Use potentially passed props
    const results = [];
    const thickness = typeof materialProps.t === 'number' ? materialProps.t : 1.0; // Needed? Not directly for stress/strain calc itself

    // Create displacement vector U (8x1 column vector)
    const U_vector = [];
    if (!displacements || displacements.length !== nodes.length) {
        console.error("Mismatch between nodes and displacements count or invalid displacements array.");
        // Fill U_vector with zeros if displacements are invalid to avoid crashing
         for (let i = 0; i < nodes.length * 2; i++) U_vector.push([0]);
    } else {
         displacements.forEach(d => {
            U_vector.push([ typeof d.x === 'number' ? d.x : 0 ]); // x displacement
            U_vector.push([ typeof d.y === 'number' ? d.y : 0 ]); // y displacement
        });
    }


    selectedGaussPoints.forEach((gaussPoint, index) => {
        const { xi, eta } = gaussPoint; // Weight not needed here

        try {
            const dNdXi = calculateDNdXi(xi, eta);
            const J = calculateJacobiMatrix(dNdXi, nodes);
            const detJ = calculateDeterminant(J);

            if (detJ <= 1e-12) {
                 results.push({
                    gaussPoint: { xi, eta },
                    error: `Invalid Jacobian determinant: ${detJ.toFixed(4)} at (xi=${xi.toFixed(4)}, eta=${eta.toFixed(4)})`
                });
                return; // Skip this Gauss point
            }

            const Jinv = calculateInverseJacobi(J, detJ);
            const dNdX = calculateDNdX(dNdXi, Jinv);
            const B = calculateBMatrix(dNdX); // 3x8

            // Calculate strain: ε = B * U (Result is 3x1)
            const strain_vector = multiplyMatrices(B, U_vector);

            // Calculate stress: σ = D * ε (Result is 3x1)
            const stress_vector = multiplyMatrices(D, strain_vector);

            // Check if matrix multiplication resulted in valid structures
            if (!strain_vector || strain_vector.length !== 3 || !stress_vector || stress_vector.length !== 3) {
                throw new Error("Matrix multiplication for strain/stress failed or yielded incorrect dimensions.");
            }

            // Extract components safely
            const strainComponents = {
                εxx: strain_vector[0]?.[0] ?? NaN,
                εyy: strain_vector[1]?.[0] ?? NaN,
                γxy: strain_vector[2]?.[0] ?? NaN  // Engineering shear strain
            };
            const stressComponents = {
                σxx: stress_vector[0]?.[0] ?? NaN,
                σyy: stress_vector[1]?.[0] ?? NaN,
                τxy: stress_vector[2]?.[0] ?? NaN
            };

            // Calculate von Mises stress
            let vonMises = NaN;
             if (Object.values(stressComponents).every(v => typeof v === 'number' && isFinite(v))) {
                 vonMises = calculateVonMisesStress(stressComponents);
             } else {
                console.warn(`[GP ${index+1}] NaN detected in stress components at (xi=${xi.toFixed(4)}, eta=${eta.toFixed(4)}). Cannot calculate von Mises stress.`);
            }


            results.push({
                gaussPoint: { xi, eta },
                strain: strainComponents,
                stress: stressComponents,
                vonMises: vonMises,
                error: null // No error for this point if we reached here
            });
        } catch(error) {
             console.error(`[GP ${index + 1}] Error during strain/stress calculation at (xi=${xi.toFixed(4)}, eta=${eta.toFixed(4)}):`, error);
             results.push({
                gaussPoint: { xi, eta },
                error: `Calculation error: ${error.message}`
            });
        }
    });

    return results;

}


export function calculateVonMisesStress(stress) {
    // Add validation for input
    if (!stress || typeof stress.σxx !== 'number' || typeof stress.σyy !== 'number' || typeof stress.τxy !== 'number') {
        return NaN;
    }
    const {σxx, σyy, τxy} = stress;
     // Check for NaN/Infinity before calculation
     if (!isFinite(σxx) || !isFinite(σyy) || !isFinite(τxy)) {
         return NaN;
     }
    // Standard formula for plane stress/strain
    const vonMisesSq = σxx*σxx - σxx*σyy + σyy*σyy + 3*τxy*τxy;
    return vonMisesSq >= 0 ? Math.sqrt(vonMisesSq) : 0; // Return 0 if result is negative due to precision issues
}

// --- Extrapolation Function (Seems OK from snippet) ---
export function extrapolateGaussToNodes(gaussResults, componentKey) {
    // ... (extrapolateGaussToNodes implementation - seems okay from snippet) ...
     // --- 1. Determine Data Path ---
    let path;
    let componentName = ''; // For logging
    if (componentKey === 'Svm') {
        path = ['vonMises'];
        componentName = 'Von Mises Stress';
    } else if (componentKey === 'Sxx') {
        path = ['stress', 'σxx'];
        componentName = 'Stress XX';
    } else if (componentKey === 'Syy') {
        path = ['stress', 'σyy'];
        componentName = 'Stress YY';
    } else if (componentKey === 'Sxy' || componentKey === 'Txy') { // Allow both common names
        path = ['stress', 'τxy'];
         componentName = 'Stress XY (Shear)';
    } else if (componentKey === 'Exx') {
        path = ['strain', 'εxx'];
        componentName = 'Strain XX';
    } else if (componentKey === 'Eyy') {
        path = ['strain', 'εyy'];
        componentName = 'Strain YY';
    } else if (componentKey === 'Exy' || componentKey === 'Gxy') { // Allow both common names
        path = ['strain', 'γxy'];
        componentName = 'Strain XY (Engineering Shear)';
    } else {
        console.error(`[Extrapolation] Unknown component key: "${componentKey}"`);
        return [0, 0, 0, 0]; // Return zero array for unknown key
    }
     console.log(`--- Extrapolating ${componentName} (${componentKey}) ---`);

    // --- 2. Initial Validation of gaussResults Structure ---
    if (!gaussResults || !Array.isArray(gaussResults)) {
        console.error(`[Extrapolation ${componentKey}] Invalid gaussResults provided (expected array). Got:`, gaussResults);
        return [0, 0, 0, 0];
    }

    const numGaussPoints = gaussResults.length;
    const defaultVal = 0; // Value to use if extraction fails or error occurs

    // --- 3. Handle based on the number of Gauss points ---

    // --- 3a. Case: 1x1 Gauss Quadrature ---
    if (numGaussPoints === 1) {
        console.log(`   Mode: 1x1 Gauss Point`);
        const result = gaussResults[0];
        let value = defaultVal;

        if (!result) {
            console.warn(`   [GP 1] Missing result object. Using default value ${defaultVal}.`);
        } else if (result.error) {
            console.warn(`   [GP 1] Calculation error reported: "${result.error}". Using default value ${defaultVal}.`);
        } else {
            // Attempt to Extract Value using Path
            let current = result;
            let validPath = true;
            for (const p of path) {
                if (current && typeof current === 'object' && Object.prototype.hasOwnProperty.call(current, p)) {
                    current = current[p];
                } else {
                    console.warn(`   [GP 1] Invalid data structure or missing path segment "${p}". Path: ${path.join('.')}. Using default value ${defaultVal}.`);
                    current = defaultVal;
                    validPath = false;
                    break;
                }
            }

            if (validPath) {
                if (typeof current === 'number' && isFinite(current)) {
                    value = current; // Success
                    console.log(`   [GP 1] Value extracted: ${value.toFixed(4)}`);
                } else {
                    console.warn(`   [GP 1] Extracted value is not a finite number (${current}). Using default value ${defaultVal}.`);
                    value = defaultVal;
                }
            }
            // else value remains defaultVal
        }

        // For 1x1, extrapolate by assigning the single value to all nodes
        const nodalValues = [value, value, value, value];
        console.log(`   Extrapolated Nodal Values: [${nodalValues.map(v=>v.toFixed(4)).join(', ')}]`);
        return nodalValues;

    // --- 3b. Case: 2x2 Gauss Quadrature ---
    } else if (numGaussPoints === 4) {
         console.log(`   Mode: 2x2 Gauss Points`);
        const gaussValues = [];
        let encounteredIssues = false;

        for (let i = 0; i < 4; i++) {
            const result = gaussResults[i];
            let value = defaultVal; // Default for this GP

            if (!result) {
                console.warn(`   [GP ${i + 1}] Missing result object. Using default value ${defaultVal}.`);
                encounteredIssues = true;
            } else if (result.error) {
                console.warn(`   [GP ${i + 1}] Calculation error reported: "${result.error}". Using default value ${defaultVal}.`);
                encounteredIssues = true;
            } else {
                // Attempt to Extract Value using Path
                let current = result;
                let validPath = true;
                for (const p of path) {
                    if (current && typeof current === 'object' && Object.prototype.hasOwnProperty.call(current, p)) {
                        current = current[p];
                    } else {
                         console.warn(`   [GP ${i + 1}] Invalid data structure or missing path segment "${p}". Path: ${path.join('.')}. Using default value ${defaultVal}.`);
                        current = defaultVal;
                        validPath = false;
                        encounteredIssues = true;
                        break;
                    }
                }

                if (validPath) {
                     if (typeof current === 'number' && isFinite(current)) {
                        value = current; // Success
                    } else {
                        console.warn(`   [GP ${i + 1}] Extracted value is not a finite number (${current}). Using default value ${defaultVal}.`);
                        value = defaultVal;
                        encounteredIssues = true;
                    }
                }
                 // else value remains defaultVal
            }
            gaussValues.push(value);
        }

        console.log(`   Extracted Gauss Point Values: [${gaussValues.map(v=>v.toFixed(4)).join(', ')}]`);
        if (encounteredIssues) {
            console.warn(`   -> Default values were used for one or more Gauss points.`);
        }

        // Perform extrapolation using the inverse shape function matrix at Gauss points
        // This specific matrix N_inv works for extrapolating 2x2 Gauss points (xi=±1/√3, eta=±1/√3)
        // to nodes (xi=±1, eta=±1) for bilinear elements.
        const N_inv = [
            [ 1.8660254037844386, -0.5                ,  0.1339745962155613 , -0.5                 ], // Node 1 (-1, -1)
            [-0.5                ,  1.8660254037844386 , -0.5                 ,  0.1339745962155613  ], // Node 2 (+1, -1)
            [ 0.1339745962155613 , -0.5                ,  1.8660254037844386 , -0.5                 ], // Node 3 (+1, +1)
            [-0.5                ,  0.1339745962155613 , -0.5                 ,  1.8660254037844386 ]  // Node 4 (-1, +1)
        ];
         // Note: Order of Gauss points in gaussResults MUST match the expected order for N_inv:
         // GP1: (-1/√3, -1/√3), GP2: (+1/√3, -1/√3), GP3: (+1/√3, +1/√3), GP4: (-1/√3, +1/√3)

        try {
            // Multiply N_inv * gaussValues_vector
             const nodalValues = multiplyMatrixVector(N_inv, gaussValues); // Assumes multiplyMatrixVector exists

             // Final check for non-finite numbers after multiplication
            if (nodalValues.some(v => typeof v !== 'number' || !isFinite(v))) {
                 console.error(`   [Extrapolation Error] Matrix multiplication resulted in non-finite values. Returning defaults. Result:`, nodalValues);
                 return [defaultVal, defaultVal, defaultVal, defaultVal];
            }
             console.log(`   Extrapolated Nodal Values: [${nodalValues.map(v=>v.toFixed(4)).join(', ')}]`);
            return nodalValues;
        } catch (e) {
            console.error(`   [Extrapolation Error] Error during matrix multiplication (N_inv * GP_values):`, e, "Input GP values:", gaussValues);
            return [defaultVal, defaultVal, defaultVal, defaultVal]; // Fallback
        }

    // --- 3c. Case: Invalid Number of Gauss Points ---
    } else {
        console.error(`[Extrapolation ${componentKey}] Invalid number of Gauss points (${numGaussPoints}). Expected 1 or 4. Got:`, gaussResults);
        return [defaultVal, defaultVal, defaultVal, defaultVal];
    }
}

// Helper for matrix-vector multiplication (needed by extrapolateGaussToNodes)
function multiplyMatrixVector(matrix, vector) {
     const numRows = matrix.length;
    if (numRows === 0) return [];
    const numCols = matrix[0].length;

    if (numCols !== vector.length) {
        throw new Error(`Matrix columns (${numCols}) must match vector length (${vector.length}) for multiplication.`);
    }

    const result = new Array(numRows).fill(0);
    for (let i = 0; i < numRows; i++) {
        for (let j = 0; j < numCols; j++) {
             const matVal = typeof matrix[i][j] === 'number' ? matrix[i][j] : 0;
             const vecVal = typeof vector[j] === 'number' ? vector[j] : 0;
             result[i] += matVal * vecVal;
        }
    }
    return result;
}