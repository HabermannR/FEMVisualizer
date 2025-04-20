import { formatter, matrixConfigs } from './matrixFormatter.js';
import {
    calculateShapeFunctions,
    calculateDNdXi,
    calculateJacobiMatrix,
    calculateDeterminant,
    calculateInverseJacobi,
    calculateDNdX,
    calculateBMatrix,
    calculateDMatrix,
    calculateStiffnessMatrix,
    reducedStiffnessMatrix,
    calculateDisplacements,
    calculateStrainAndStress,
    calculateVonMisesStress
} from './calculation.js';

// Gauss point configurations
const gaussPoints1D = {
    1: {
        points: [[0.0]],
        weights: [2.0]
    },
    2: {
        points: [[-0.577350269189626], [0.577350269189626]],
        weights: [1.0, 1.0]
    },
    3: {
        points: [[-0.774596669241483], [0.0], [0.774596669241483]],
        weights: [0.555555555555556, 0.888888888888889, 0.555555555555556]
    }
};

const gaussPoints2D = {
    1: {
        points: [[0.0, 0.0]],
        weights: [4.0]
    },
    2: {
        points: [
            [-0.577350269189626, -0.577350269189626],
            [-0.577350269189626, 0.577350269189626],
            [0.577350269189626, -0.577350269189626],
            [0.577350269189626, 0.577350269189626]
        ],
        weights: [1.0, 1.0, 1.0, 1.0]
    }
};

function updateSection(sectionId, content) {
    const section = document.getElementById(sectionId);
    if (section) {
        section.innerHTML = content;
    } else {
        console.warn(`Section with id ${sectionId} not found`);
    }
}

export function updateLessons(currentGaussPoint, nodes, forceArrow, activeQuadratureRuleName) {
    const {xi, eta} = currentGaussPoint;
    
    const N = calculateShapeFunctions(xi, eta);
    const dNdXi = calculateDNdXi(xi, eta);
    const X = nodes.map(node => [node.x, node.y]);
    const J = calculateJacobiMatrix(dNdXi, nodes);
    const detJ = calculateDeterminant(J);
    const Jinv = calculateInverseJacobi(J, detJ);
    const dNdX = calculateDNdX(dNdXi, Jinv);
    const B = calculateBMatrix(dNdX);
    const D = calculateDMatrix();
    const stiffnessData = calculateStiffnessMatrix(nodes, undefined, activeQuadratureRuleName);
    const reducedstiffnessData = reducedStiffnessMatrix(nodes, undefined, activeQuadratureRuleName);
    const displacements = calculateDisplacements(nodes, forceArrow, undefined, activeQuadratureRuleName);

    // Symbolic matrices definitions
    const symbolicN = [
        ['¼(1-ξ)(1-η)', '¼(1+ξ)(1-η)', '¼(1+ξ)(1+η)', '¼(1-ξ)(1+η)']
    ];

    const symbolicDNdXi = [
        ['∂N₁/∂ξ', '∂N₁/∂η'],
        ['∂N₂/∂ξ', '∂N₂/∂η'],
        ['∂N₃/∂ξ', '∂N₃/∂η'],
        ['∂N₄/∂ξ', '∂N₄/∂η']
    ];

    const symbolicDNdXiRes = [
        ['-¼(1-η)', '-¼(1-ξ)'],
        ['¼(1-η)',  '-¼(1+ξ)'],
        ['¼(1+η)',  '¼(1+ξ)'],
        ['-¼(1+η)', '¼(1-ξ)']
    ];

    const symbolicJ = [
        ['∂x/∂ξ', '∂y/∂ξ'],
        ['∂x/∂η', '∂y/∂η']
    ];

    const symbolicP = [
        ['∂N₁/∂ξ', '∂N₂/∂ξ', '∂N₃/∂ξ', '∂N₄/∂ξ'],
        ['∂N₁/∂η', '∂N₂/∂η', '∂N₃/∂η', '∂N₄/∂η']
    ];

    const symbolicX = [
        ['x₁', 'y₁'],
        ['x₂', 'y₂'],
        ['x₃', 'y₃'],
        ['x₄', 'y₄']
    ];

    const symbolicJinv = [
        [' ∂y/∂η/|J|', '-∂y/∂ξ/|J|'],
        ['-∂x/∂η|J|', ' ∂x/∂ξ/|J|']
    ];

    const symbolicDNdX = [
        ['∂N₁/∂x', '∂N₁/∂y'],
        ['∂N₂/∂x', '∂N₂/∂y'],
        ['∂N₃/∂x', '∂N₃/∂y'],
        ['∂N₄/∂x', '∂N₄/∂y']
    ];

    const symbolicB = [
        ['∂N₁/∂x', '0', '∂N₂/∂x', '0', '∂N₃/∂x', '0', '∂N₄/∂x', '0'],
        ['0', '∂N₁/∂y', '0', '∂N₂/∂y', '0', '∂N₃/∂y', '0', '∂N₄/∂y'],
        ['∂N₁/∂y', '∂N₁/∂x', '∂N₂/∂y', '∂N₂/∂x', '∂N₃/∂y', '∂N₃/∂x', '∂N₄/∂y', '∂N₄/∂x']
    ];
    const symbolicD = [
        ['1', 'ν', '0'],
        ['ν', '1', ' 0'],
        ['0', '0', '(1-ν)/2']
    ];


    if (!document.getElementById('fem-calculations')) {
        document.getElementById('shapeFunctions').innerHTML = `
            <h3>Finite Element Analysis Step-by-Step Calculation</h3>
            <div id="fem-calculations">
                <details id="intro-section" open>
                    <summary><h4>1. Introduction and Node Coordinates</h4></summary>
                    <div id="intro-content" class="section-content"></div>
                </details>

                <details id="gauss-section">
                    <summary><h4>2. Gauss Point Location</h4></summary>
                    <div id="gauss-content" class="section-content"></div>
                </details>

                <details id="shape-section">
                    <summary><h4>3. Shape Functions</h4></summary>
                    <div id="shape-content" class="section-content"></div>
                </details>

                <details id="derivatives-section">
                    <summary><h4>4. Shape Function Derivatives</h4></summary>
                    <div id="derivatives-content" class="section-content"></div>
                </details>

                <details id="jacobian-section">
                    <summary><h4>5. Jacobian Matrix</h4></summary>
                    <div id="jacobian-content" class="section-content"></div>
                </details>

                <details id="inverse-jacobian-section">
                    <summary><h4>6. Inverse Jacobian</h4></summary>
                    <div id="inverse-jacobian-content" class="section-content"></div>
                </details>

                <details id="physical-derivatives-section">
                    <summary><h4>7. Physical Derivatives</h4></summary>
                    <div id="physical-derivatives-content" class="section-content"></div>
                </details>

                <details id="b-matrix-section">
                    <summary><h4>8. B-Matrix</h4></summary>
                    <div id="b-matrix-content" class="section-content"></div>
                </details>

                <details id="stiffness-section">
                    <summary><h4>9. Stiffness Matrix</h4></summary>
                    <div id="stiffness-content" class="section-content"></div>
                </details>

                <details id="reduced-section">
                    <summary><h4>10. Stiffness Matrix reduction</h4></summary>
                    <div id="reduced-content" class="section-content"></div>
                </details>

                <details id="disp-section">
                    <summary><h4>11. Solving for the node displacements</h4></summary>
                    <div id="disp-content" class="section-content"></div>
                </details>

                <details id="strain-stress-section">
                    <summary><h4>12. Strain and Stress Analysis</h4></summary>
                    <div id="strain-stress-content" class="section-content"></div>
                </details>
            </div>
        `;
    }

    // Update each section's content
    updateSection('intro-content', `
        <p>The finite element method (FEM) discretizes complex geometries into simpler elements for numerical analysis.
        For a 4-node quadrilateral element:</p>
        ${formatter.createMatrix(
            nodes.map(node => [node.x, node.y]),
            {
                colHeaders: ['X', 'Y'],
                rowHeaders: nodes.map((_, i) => `Node ${i + 1}`),
                title: 'Node Coordinates'
            }
        )}
        <p>You can move the nodes and the force arrow in the left Physical Element view, and the middle Gauss point in the middle
        <p>Attention, for the graphics, the y axis is inverted, meaning -1, -1 is the left, top corner!</p>
    `);

    updateSection('gauss-content', `
        <p>Gauss-Legendre quadrature is a numerical integration method characterized by:</p>
        <ul>
            <li>Strategically positioned sampling points that optimize integration accuracy</li>
            <li>Associated weight factors for each sampling point</li>
            <li>Natural coordinates defined in the range [-1, +1]</li>
        </ul>

        <h5>1D Integration Points</h5>
        <div class="matrix-equation-container">
            ${Object.entries(gaussPoints1D).map(([order, data]) =>
                formatter.createMatrix(
                    data.points.map((point, idx) => [
                        point[0],
                        data.weights[idx]
                    ]),
                    {
                        title: `${order}-Point Rule`,
                        colHeaders: ['Position (ξ)', 'Weight'],
                        rowHeaders: data.points.map((_, i) => `Point ${i + 1}`),
                        customClasses: ['gauss-table']
                    }
                )
            ).join('')}
        </div>
        <h5>2D Integration Points</h5>
        <div class="matrix-equation-container">
            ${Object.entries(gaussPoints2D).map(([order, data]) =>
                formatter.createMatrix(
                    data.points.map((point, idx) => [
                        point[0],
                        point[1],
                        data.weights[idx]
                    ]),
                    {
                        title: `${order}-Point Rule`,
                        colHeaders: ['ξ', 'η', 'Weight'],
                        rowHeaders: data.points.map((_, i) => `Point ${i + 1}`),
                        customClasses: ['gauss-table']
                    }
                )
            ).join('')}
        </div>
        <p><em>Note: For 2D, the points form a grid pattern, resulting in 4 points (2×2 grid).</em></p>

        <h5>Current Integration Point</h5>
        ${formatter.createMatrix([[xi, eta]], {
            colHeaders: ['ξ', 'η'],
            customClasses: ['current-point']
        })}
    `);

    updateSection('shape-content', `
        <p>Shape functions are fundamental to FEM, providing interpolation within the element:</p>
        <ul>
            <li>They determine how field variables vary within the element</li>
            <li>Each shape function equals 1 at its associated node and 0 at other nodes</li>
            <li>The sum of all shape functions equals 1 at any point (partition of unity)</li>
        </ul>

        ${formatter.createMatrix(symbolicN, {
            ...matrixConfigs.N,
            title: 'Analytical Form'
        })}
        <p></p>
        ${formatter.createMatrix([N], {
            ...matrixConfigs.N,
            title: `Numerical values at (ξ=${xi.toFixed(2)}, η=${eta.toFixed(2)})`
        })}
    `);

    updateSection('derivatives-content', `
        <p>Derivatives of shape functions with respect to natural coordinates (ξ,η) are needed to
        compute the Jacobian matrix and ultimately relate natural to physical coordinates.</p>

        <div class="matrix-equation-container">
            ${formatter.createMatrix(symbolicDNdXi, {
                ...matrixConfigs['dN/dξ'],
                title: 'Analytical Form'
            })}

            ${formatter.createMatrix(symbolicDNdXiRes, {
                ...matrixConfigs['dN/dξ'],
                title: 'Analytical Form resolved'
            })}

            ${formatter.createMatrix(dNdXi, {
                ...matrixConfigs['dN/dξ'],
                title: 'Numerical Values'
            })}
        </div>
    `);

    updateSection('jacobian-content', `
        <p>The Jacobian matrix (J) maps derivatives from natural (ξ,η) to physical (x,y) coordinates.
        It represents the local linear transformation between coordinate systems.</p>

        <p>J = P × X, where:</p>
        <div class="matrix-equation-container">
            ${formatter.createMatrix(symbolicJ, {
                ...matrixConfigs.J,
                title: 'Analytical Form (J)'
            })} =
            ${formatter.createMatrix(symbolicP, {
                title: 'P (Shape Function Derivatives)'
            })} ×
            ${formatter.createMatrix(symbolicX, {
                title: 'X (Nodal Coordinates)'
            })}
        </div>

        <p>Numerical values:</p>
        <div class="matrix-equation-container">
            ${formatter.createMatrix(J, {
                ...matrixConfigs.J,
                title: 'J'
            })} =
            ${formatter.createMatrix(dNdXi, {
                ...matrixConfigs['dN/dξ'],
                title: 'P (Shape Function Derivatives)'
            })} ×
            ${formatter.createMatrix(
                nodes.map(node => [node.x, node.y]),
                {
                    colHeaders: ['X', 'Y'],
                    rowHeaders: nodes.map((_, i) => `Node ${i + 1}`),
                    title: 'Node Coordinates'
                }
            )}
        </div>
    `);

     updateSection('inverse-jacobian-content', `
        <p>The inverse Jacobian (J⁻¹) is used to transform derivatives from natural to physical
        coordinates. It's essential for computing strains in physical space.</p>

        ${formatter.createMatrix(symbolicJinv, {
            ...matrixConfigs['J⁻¹'],
            title: 'Analytical Form'
        })}

        <p>With determinant |J| = ${detJ.toFixed(6)}</p>

        ${formatter.createMatrix(Jinv, {
            ...matrixConfigs['J⁻¹'],
            title: 'Numerical Values'
        })}
    `);

    updateSection('physical-derivatives-content', `
        <p>The derivatives with respect to physical coordinates (x, y) are obtained from the derivatives with respect to natural coordinates (ξ, η) using the inverse Jacobian matrix J⁻¹. For each shape function Nᵢ:</p>
        <div class="matrix-equation-container">
            ${formatter.createMatrix(symbolicDNdX, {
                ...matrixConfigs['dN/dx'],
                title: 'Analytical Form dN/dx'
            })} = J⁻¹ x

            ${formatter.createMatrix(symbolicDNdXi, {
                ...matrixConfigs['dN/dξ'],
                title: 'P (Shape Function Derivatives)'
            })}
        </div>
        <p></p>
        <div class="matrix-equation-container">
            ${formatter.createMatrix(dNdX, {
                ...matrixConfigs['dN/dx'],
                title: 'Numerical Values'
            })} =
            ${formatter.createMatrix(Jinv, {
                ...matrixConfigs['J⁻¹'],
                title: 'Numerical Values'
            })}
            x

            ${formatter.createMatrix(dNdXi, {
                ...matrixConfigs['dN/dξ'],
                title: 'P (Shape Function Derivatives)'
            })}
        </div>
    `);


    updateSection('b-matrix-content', `
        <p>The B-matrix relates nodal displacements to strains. It's constructed from the physical
        derivatives and is crucial for computing the element stiffness matrix. It is also called
        the strain-displacement matrix.</p>

        <div class="formula">
            {ε} = [B]{u}
        </div>

        <p>where {ε} = [εxx εyy γxy]ᵀ and {u} contains nodal displacements</p>

        ${formatter.createMatrix(symbolicB, {
            ...matrixConfigs.B,
            title: 'Analytical Form'
        })}

        ${formatter.createMatrix(B, {
            ...matrixConfigs.B,
            title: 'Numerical Values'
        })}
    `);

    updateSection('stiffness-content', `
        <p>The <strong>Element Stiffness Matrix (K)</strong> is a fundamental concept in FEM. It establishes the
        linear relationship between the forces applied at the element's nodes and the resulting displacements
        of those nodes: <strong>{F} = [K]{u}</strong>.</p>

        <p>It is computed by integrating the product of the transposed Strain-Displacement matrix (B<sup>T</sup>),
        the Material Constitutive matrix (D), and the Strain-Displacement matrix (B) over the element's area (or volume in 3D).
        This integration effectively sums up the stiffness contributions from every infinitesimal point within the element.</p>

        <div class="formula">
             K = ∫<sub>A</sub> B(x, y)<sup>T</sup> D B(x, y) dA
        </div>

        <p>However, calculations are typically performed in the simpler <strong>Natural Coordinate System (ξ, η)</strong> ranging from -1 to +1.
        The integral transforms to:</p>

        <div class="formula">
             K = ∫<sub>-1</sub><sup>1</sup> ∫<sub>-1</sub><sup>1</sup> B(ξ, η)<sup>T</sup> D B(ξ, η) |J(ξ, η)| dξ dη * t
        </div>
         <p>Where 't' is the element thickness (assumed = 1.0 if not specified).</p>

        <p>Where the components are:</p>
        <ul>
            <li><strong>B(ξ, η)</strong>: The Strain-Displacement Matrix (calculated in the previous step). It links nodal displacements {u} to the strain {ε} at any point (ξ, η) within the element. Notice it depends on the location (ξ, η).</li>
            <li><strong>D</strong>: The Material Constitutive Matrix. It relates stress {σ} to strain {ε} ({σ} = [D]{ε}) based on material properties (Young's Modulus E, Poisson's Ratio ν). For this element, we assume linear elastic behavior and plane stress conditions.</li>
            <li><strong>B(ξ, η)<sup>T</sup></strong>: The transpose of the B matrix. It is used (conceptually via virtual work) to transform stresses back into equivalent nodal forces.</li>
            <li><strong>|J(ξ, η)|</strong>: The Determinant of the Jacobian Matrix (calculated earlier). It is the scaling factor needed to correctly relate an infinitesimal area <code>dξ dη</code> in the natural coordinate system to the corresponding area <code>dx dy</code> in the physical system. It accounts for the element's shape and distortion. This also generally depends on (ξ, η).</li>
            <li><strong>t</strong>: Element thickness.</li>
        </ul>

        <h5>Numerical Integration: Gaussian Quadrature</h5>
        <p>Because the B matrix and the Jacobian determinant |J| typically vary across the element, this integral is usually evaluated numerically using <strong>Gaussian Quadrature</strong>. This method approximates the integral as a weighted sum of the integrand evaluated at specific points called <strong>Gauss Points</strong> within the natural coordinate domain [-1, 1] x [-1, 1].</p>
        <p>The integral is approximated as:</p>
        <div class="formula">
             K ≈ Σ<sub>i=1</sub><sup>N<sub>gp</sub></sup> [ B(ξ<sub>i</sub>, η<sub>i</sub>)<sup>T</sup> D B(ξ<sub>i</sub>, η<sub>i</sub>) |J(ξ<sub>i</sub>, η<sub>i</sub>)| * w<sub>i</sub> * t ]
        </div>
        <p>Where:</p>
        <ul>
            <li>N<sub>gp</sub> is the number of Gauss points (1 for '1x1' rule, 4 for '2x2' rule).</li>
            <li>(ξ<sub>i</sub>, η<sub>i</sub>) are the coordinates of the i-th Gauss point.</li>
            <li>w<sub>i</sub> is the weight associated with the i-th Gauss point.</li>
            <li>The term inside the square brackets <code>[...]</code> is the <strong>contribution</strong> of the i-th Gauss point to the total stiffness matrix.</li>
        </ul>
        <p>Currently using the <strong>${activeQuadratureRuleName}</strong> Gauss quadrature rule.</p>

        <h5>Material Matrix (D)</h5>
        <p>The Material Matrix D for plane stress, assuming isotropic material, is given by:</p>
        <div class="formula"> D = (E / (1 - ν²)) * [ 1   ν   0 ; ν   1   0 ; 0   0  (1-ν)/2 ]</div>

        <div class="matrix-equation-container">
            ${formatter.createMatrix(symbolicD, {
                ...matrixConfigs.D,
                title: 'Symbolic Form (Factor E/(1-ν²)/t needs to be applied)' // Added thickness
            })}

            ${formatter.createMatrix(D, {
                ...matrixConfigs.D,
                title: 'Numerical Values (Using E=210000 MPa, ν=0.3)'
            })}
        </div>

        <h5>Calculating the Element Stiffness Matrix</h5>

        ${(() => {
            // --- Show Partial Contributions based on rule ---
            if (activeQuadratureRuleName === '2x2' && stiffnessData.contributions && stiffnessData.contributions.length === 4) {
                return `
                    <h6>Partial Stiffness Matrices (2x2 Rule Contributions)</h6>
                    <p>Below are the calculated contributions from each of the 4 Gauss points. Each contribution is an 8x8 matrix calculated as: <code>K<sub>gp, i</sub> = B<sub>i</sub><sup>T</sup> D B<sub>i</sub> |J<sub>i</sub>| w<sub>i</sub> t</code> (w<sub>i</sub>=1 for all points in 2x2 rule)</p>
                    <div class="gauss-contributions">
                        ${stiffnessData.contributions.map((contrib, i) => `
                            <div class="gauss-contribution">
                                <h6>Contribution from Gauss Point ${i + 1} (ξ=${contrib.gaussPoint.xi.toFixed(4)}, η=${contrib.gaussPoint.eta.toFixed(4)}, w=1)</h6>
                                <p>Determinant |J| at this point: ${contrib.detJ.toFixed(4)}</p>
                                ${contrib.error ? `<p class="error-message">Error: ${contrib.error}</p>` : ''}
                                ${formatter.createMatrix(contrib.contribution || Array(8).fill(0).map(() => Array(8).fill(0)), { // Provide default if contribution missing
                                    colHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
                                    rowHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
                                    isStiffness: true,
                                    title: `K<sub>Gauss ${i + 1}</sub>`,
                                    customClasses: ['contribution-matrix']
                                })}
                            </div>
                        `).join('')}
                    </div>
                `;
            } else if (activeQuadratureRuleName === '1x1' && stiffnessData.contributions && stiffnessData.contributions.length === 1) {
                 const contrib = stiffnessData.contributions[0];
                 return `
                     <h6>Stiffness Matrix Contribution (1x1 Rule)</h6>
                     <p>With the 1x1 rule, the stiffness is evaluated only at the center (ξ=0, η=0) with weight w=4: <code>K<sub>gp</sub> = B<sub>0</sub><sup>T</sup> D B<sub>0</sub> |J<sub>0</sub>| * 4 * t</code></p>
                     <div class="gauss-contributions">
                         <div class="gauss-contribution">
                              <h6>Contribution from Gauss Point 1 (ξ=0, η=0, w=4)</h6>
                              <p>Determinant |J| at center: ${contrib.detJ.toFixed(4)}</p>
                              ${contrib.error ? `<p class="error-message">Error: ${contrib.error}</p>` : ''}
                              ${formatter.createMatrix(contrib.contribution || Array(8).fill(0).map(() => Array(8).fill(0)), {
                                    colHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
                                    rowHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
                                    isStiffness: true,
                                    title: `K<sub>Gauss 1</sub>`,
                                    customClasses: ['contribution-matrix']
                              })}
                         </div>
                     </div>
                `;
            } else {
                // Handle cases where contributions might be missing or length doesn't match rule
                 console.warn(`Stiffness contributions array mismatch or missing. Rule: ${activeQuadratureRuleName}, Contributions:`, stiffnessData.contributions);
                 return '<p>Stiffness contributions data is unavailable or inconsistent with the selected rule.</p>';
            }
        })()}

        <!-- Hourglass Control Explanation -->
        <h6>Hourglass Control (for 1x1 Reduced Integration)</h6>
        <p>When using <strong>reduced integration</strong> schemes like the 1x1 rule for a 4-node quadrilateral, the single integration point may not be sufficient to accurately capture all possible deformation modes of the element. Specifically, certain deformation patterns, known as <strong>hourglass modes</strong> or <strong>zero-energy modes</strong>, can occur without causing any strain at the element's center (the single Gauss point). This means the standard stiffness matrix calculation using the 1x1 rule will assign zero stiffness (zero energy) to these modes, potentially leading to unrealistic, jagged, or oscillatory deformation patterns in the overall structure, especially in mesh regions dominated by bending.</p>
        <p>To counteract this, an <strong>hourglass stabilization</strong> technique is often employed. This involves adding a small, artificial stiffness specifically designed to penalize these hourglass modes. A common approach involves:</p>
        <ul>
            <li>Defining an <strong>hourglass vector</strong> (γ) that characterizes the hourglass shape (e.g., [1, -1, 1, -1]).</li>
            <li>Projecting out the rigid body and constant strain components from this vector to get a <strong>corrected gamma vector</strong> (γ̂) that is orthogonal to these modes.</li>
            <li>Constructing an <strong>hourglass stiffness matrix (K<sub>hg</sub>)</strong>, typically proportional to the element's shear modulus (G), area (A), thickness (t), and a user-defined stabilization parameter (α):
                <div class="formula"> K<sub>hg</sub> ≈ α * G * A * t * (γ̂ ⊗ γ̂) ⊗ I₂ </div>
                (Where ⊗ is an outer product and I₂ is the 2x2 identity matrix).
            </li>
            <li>Adding this K<sub>hg</sub> to the stiffness matrix calculated at the single Gauss point: <strong>K<sub>total</sub> = K<sub>Gauss 1</sub> + K<sub>hg</sub></strong>.</li>
        </ul>
        <p>The stabilization parameter <strong>α (alpha)</strong> controls the amount of artificial stiffness added. It needs careful selection – too small, and hourglassing might still occur; too large, and the element becomes overly stiff ("locking"). A typical range might be 0.01 to 0.1.</p>

        ${(() => {
            // --- Show Hourglass Details if 1x1 rule was used and stabilization applied ---
            if (activeQuadratureRuleName === '1x1' && stiffnessData.hgDetails) {
                const details = stiffnessData.hgDetails;
                //removed materialprops.t
                const alpha = stiffnessData.hgDetails?.hgStiffnessMagnitude && stiffnessData.hgDetails?.G && stiffnessData.hgDetails?.Area
                              ? (stiffnessData.hgDetails.hgStiffnessMagnitude / (stiffnessData.hgDetails.G * stiffnessData.hgDetails.Area * (1.0))).toFixed(3)
                              : "N/A (Calculation Error?)"; // Calculate alpha back if possible

                let hgContent = `
                    <div class="hourglass-details">
                        <h6>Applied Hourglass Stabilization Details (α ≈ ${alpha}):</h6>
                `;
                if (details.error) {
                    hgContent += `<p class="error-message"><strong>Warning/Error during stabilization:</strong> ${details.error}</p>`;
                }
                if (details.K_hg) {
                     hgContent += `
                        <p>The following Hourglass Stiffness Matrix (K<sub>hg</sub>) was calculated and added:</p>
                        ${formatter.createMatrix(details.K_hg, {
                            colHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
                            rowHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
                            isStiffness: true,
                            title: `K<sub>hg</sub>`,
                            customClasses: ['contribution-matrix', 'hourglass-matrix']
                        })}
                     `;
                } else {
                    hgContent += `<p>Hourglass stiffness matrix (K<sub>hg</sub>) was not calculated or added, likely due to an error (e.g., distorted element at center).</p>`;
                }
                hgContent += `</div>`;
                return hgContent;
            } else if (activeQuadratureRuleName === '1x1') {
                 return `<p>Hourglass stabilization was not applied (alpha might be zero or calculation skipped).</p>`;
            } else {
                 // Only mention it applies to 1x1 if currently in 2x2 mode
                return `<p>Hourglass control is typically relevant for the 1x1 reduced integration rule, not the currently selected 2x2 rule.</p>`;
            }
        })()}


        <h6>Complete Element Stiffness Matrix (K<sub>total</sub>)</h6>
        <p>The final Element Stiffness Matrix K<sub>total</sub> is obtained by summing the contributions from all Gauss points ${activeQuadratureRuleName === '1x1' && stiffnessData.hgDetails && !stiffnessData.hgDetails.error ? ' and adding the hourglass stabilization matrix K<sub>hg</sub>' : activeQuadratureRuleName === '1x1' ? '' : '(K<sub>total</sub> = Σ K<sub>Gauss i</sub>)'}.</p>

        ${formatter.createMatrix(stiffnessData.K, { // K already includes K_hg if calculated
            colHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
            rowHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
            title: 'Complete Stiffness Matrix (K<sub>total</sub>)',
            isStiffness: true,
            customClasses: ['stiffness-matrix']
        })}

        <div class="matrix-explanation">
            <h6>Matrix Structure Reminder:</h6>
            <p>Each row and column corresponds to a specific nodal Degree of Freedom (DOF):</p>
            <ul>
                <li>Rows/Columns 1-2: Node 1 displacements (u₁, v₁)</li>
                <li>Rows/Columns 3-4: Node 2 displacements (u₂, v₂)</li>
                <li>Rows/Columns 5-6: Node 3 displacements (u₃, v₃)</li>
                <li>Rows/Columns 7-8: Node 4 displacements (u₄, v₄)</li>
            </ul>
            <p>The value K<sub>ij</sub> represents the force required at DOF <i>i</i> to produce a unit displacement at DOF <i>j</i> (while all other DOFs are held fixed).</p>
        </div>
    `);

    updateSection('reduced-content', `
        <p>The boundary conditions, fixed in x and y direction on node 4 (DOFs 6, 7) and fixed in y direction on node 3 (DOF 5), reduce the degrees of freedom considered in the global system. The rows and columns corresponding to these fixed DOFs are removed from the global stiffness matrix and force vector.</p>
        <h5>Reduced Element Stiffness Matrix (K<sub>r</sub>)</h5>
        <p>This shows the 5x5 submatrix corresponding to the free DOFs (1, 2, 3, 4, 5):</p>
        ${formatter.createMatrix(reducedstiffnessData, {
            colHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃' ],
            rowHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃'],
            title: 'Reduced Stiffness Matrix (K<sub>r</sub>)',
            isStiffness: true,
            customClasses: ['stiffness-matrix']
        })}
    `);

     // --- Start of Added Code ---
    // Prepare data for displacement field calculation display
    const dVector = displacements.reduce((acc, disp) => {
        acc.push([disp.x]);
        acc.push([disp.y]);
        return acc;
    }, []); // Flatten nodal displacements into [u1, v1, u2, v2, ...] column vector

    const N_matrix = [
        [N[0], 0, N[1], 0, N[2], 0, N[3], 0],
        [0, N[0], 0, N[1], 0, N[2], 0, N[3]]
    ]; // Expanded shape function matrix [2x8]

    // Calculate displacement field {u} at the current Gauss point (xi, eta)
    const displacementField = [ [0], [0] ]; // Initialize result [ux, uy] column vector
    for (let i = 0; i < 8; i++) {
        displacementField[0][0] += N_matrix[0][i] * dVector[i][0]; // ux = N1*u1 + N2*u2 + ...
        displacementField[1][0] += N_matrix[1][i] * dVector[i][0]; // uy = N1*v1 + N2*v2 + ...
    }

     updateSection('disp-content', `
        <p>After applying boundary conditions (fixing DOFs 5, 6, 7) and solving the reduced system K<sub>r</sub> * u<sub>r</sub> = F<sub>r</sub>, we obtain the displacements for the free degrees of freedom (u<sub>r</sub>). These are then assembled back into the full displacement vector {u}, setting fixed DOFs to zero.</p>
        <div class="matrix-equation-container">
             ${formatter.createMatrix(
                 nodes.map(node => [node.x, node.y]),
                 { /* config */ }
             )}
             ${formatter.createMatrix(
                 displacements.map(d => [d.x, d.y]),
                 { /* config */ }
             )}
            ${formatter.createMatrix(
                 nodes.map((node, i) => [node.x + displacements[i].x, node.y + displacements[i].y]),
                 { /* config */ }
            )}
        </div>
        <p>Note:</p>
        <ul>
            <li>Node 4 is fixed (u₄=0, v₄=0) -> DOFs 6, 7</li>
            <li>Node 3 can only move horizontally (v₃=0, roller support) -> DOF 5</li>
            <li>Force is applied at Node 2 with magnitude ${Math.round(forceArrow.length)}N at angle ${Math.round(forceArrow.angle * 180 / Math.PI)}° (Affects DOFs 2, 3)</li>
        </ul>
        <hr>
        <h4>Calculating Displacement Field within the Element</h4>
        <p>While the previous section shows the calculated displacements at the nodes, the shape functions allow us to interpolate the displacement field {u} = [uₓ, u\u1D67]ᵀ at any point (ξ, η) inside the element using the nodal displacements {d}.</p>
        <div class="formula">{u(ξ,η)} = [N<sub>matrix</sub>(ξ,η)] {d}</div>
        <p>Where {d} is the vector of nodal displacements: {d} = [u₁, v₁, u₂, v₂, u₃, v₃, u₄, v₄]ᵀ, and [N<sub>matrix</sub>] is constructed from the shape functions:</p>
        <div class="matrix-equation-container">
             ${formatter.createMatrix([ /* Symbolic N_matrix */ ], { /* config */ })}
        </div>
        <p>Using the calculated nodal displacements and evaluating the shape functions at the current Gauss point (ξ=${xi.toFixed(3)}, η=${eta.toFixed(3)}):</p>
        <div class="matrix-equation-container">
             ${formatter.createMatrix(displacementField, { /* config */ })} =
             ${formatter.createMatrix(N_matrix, { /* config */ })} ×
             ${formatter.createMatrix(dVector, { /* config */ })}
        </div>
        <p>This calculation shows how the displacement varies *within* the element based on the nodal values. This interpolated displacement field is then used to calculate strains and stresses (as shown in the next section, where the strain calculation uses {ε} = [B]{d}).</p>
    `);

    // --- Strain/Stress Section Update ---
    updateSection('strain-stress-content', `
        <p>Using the nodal displacement results {d}, we can calculate strains and stresses at each Gauss point used in the <strong>${activeQuadratureRuleName}</strong> integration rule:</p>

        <div class="formula">
            Strain at GP<sub>i</sub>: {ε}<sub>i</sub> = [B(ξ<sub>i</sub>, η<sub>i</sub>)]{d} <br>
            Stress at GP<sub>i</sub>: {σ}<sub>i</sub> = [D]{ε}<sub>i</sub>
        </div>

        <div class="formula-explanation">
            <p>Where:</p>
            <ul>
                <li>{ε} = [εxx εyy γxy]ᵀ (normal and engineering shear strains)</li>
                <li>{σ} = [σxx σyy τxy]ᵀ (normal and shear stresses)</li>
                <li>[B(ξ<sub>i</sub>, η<sub>i</sub>)] is the B-Matrix evaluated at the specific Gauss point.</li>
                 <li>{d} is the full vector of calculated nodal displacements [u₁, v₁, u₂, v₂, u₃, v₃, u₄, v₄]ᵀ.</li>
            </ul>
        </div>

        ${(() => {
            // Calculate using the active rule name
            const results = calculateStrainAndStress(nodes, displacements, undefined, activeQuadratureRuleName);
            if (!results || results.length === 0) {
                return '<p>No strain/stress results calculated.</p>';
            }

            // --- Display results for each Gauss point (1 or 4) ---
            return results.map((result, idx) => {
                // Check for explicit error property
                if (result.error) {
                    return `
                        <div class="gauss-point-results gauss-point-error">
                            <h5>Gauss Point ${idx + 1} (ξ=${result.gaussPoint?.xi?.toFixed(4) ?? 'N/A'}, η=${result.gaussPoint?.eta?.toFixed(4) ?? 'N/A'})</h5>
                            <div class="error-message">
                                Calculation Error: ${result.error}
                            </div>
                        </div>
                    `;
                }
                // Check for expected data structure if no error reported
                else if (result.gaussPoint && result.strain && result.stress && typeof result.vonMises === 'number') {
                    // Format valid results
                    return `
                        <div class="gauss-point-results">
                            <h5>Gauss Point ${idx + 1} (ξ=${result.gaussPoint.xi.toFixed(4)}, η=${result.gaussPoint.eta.toFixed(4)})</h5>
                            <div class="matrix-equation-container">
                                ${formatter.createMatrix([
                                    [result.strain.εxx],
                                    [result.strain.εyy],
                                    [result.strain.γxy]
                                ], {
                                    rowHeaders: ['εxx', 'εyy', 'γxy'], // Use gamma for shear strain
                                    title: 'Strain Components',
                                    precision: 6,
                                    customClasses: ['matrix-container']
                                })}

                                ${formatter.createMatrix([
                                    [result.stress.σxx],
                                    [result.stress.σyy],
                                    [result.stress.τxy]
                                ], {
                                    rowHeaders: ['σxx', 'σyy', 'τxy'], // Use tau for shear stress
                                    title: 'Stress Components (MPa)',
                                    precision: 1,
                                    customClasses: ['matrix-container']
                                })}
                            </div>

                            <div class="von-mises-value">
                                von Mises Stress: ${result.vonMises.toFixed(2)} MPa
                            </div>
                        </div>
                    `;
                }
                // Handle cases where structure is unexpected but no error was explicitly set
                else {
                     return `
                        <div class="gauss-point-results gauss-point-warning">
                            <h5>Gauss Point ${idx + 1} (ξ=${result.gaussPoint?.xi?.toFixed(4) ?? 'N/A'}, η=${result.gaussPoint?.eta?.toFixed(4) ?? 'N/A'})</h5>
                            <div class="warning-message">
                                Warning: Incomplete or unexpected data format for this point.
                            </div>
                        </div>
                    `;
                }
            }).join('');
        })()}

        <div class="stress-explanation">
            <h5>Understanding the Results:</h5>
            <ul>
                <li>εxx, εyy: Normal strains in x and y directions</li>
                <li>γxy: Engineering shear strain (change in angle)</li>
                <li>σxx, σyy: Normal stresses (force/area) in x and y directions</li>
                <li>τxy: Shear stress</li>
                <li>von Mises stress: An equivalent stress measure used to predict yielding of ductile materials under complex loading. It combines normal and shear stresses into a single scalar value.</li>
            </ul>
            <p>Note: The values shown are calculated at the specific Gauss points. The contour plot in the 'Result' view extrapolates these values to the nodes and interpolates across the element.</p>
        </div>
    `);
}