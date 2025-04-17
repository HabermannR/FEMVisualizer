import { updateLessons } from './updateLessons.js';
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
    calculateDisplacements,
    calculateStrainAndStress,
    extrapolateGaussToNodes,
    standardGaussPoints
} from './calculation.js';


// reenable element distortion metric
// todo Abaqus export?


// Constants
const PHYSICAL_SPACE = {
    width: 1,
    height: 1
};

const DISPLACEMENT_SCALE_FACTOR = 30;

const COLORS = {
    element: '#2196F3',    // Material Blue
    node: '#F44336',       // Material Red
    text: '#212121',       // Dark Gray
    debug: '#757575',      // Medium Gray
    boundaryCondition: '#424242',    // Material Dark Gray
    boundaryConditionFill: '#FFFFFF', // Material White
    force: '#4CAF50',      // Material Green
    forceSelected: '#81C784' // Material Light Green
};

const NODE_RADIUS = 8;
const NODE_OUTLINE_COLOR = 'black';
const NODE_OUTLINE_WIDTH = 2;
const NODE_LABEL_FONT = '14px Arial';
const NODE_LABEL_COLOR = 'black';

const MARGIN = 40;

// Initial node positions in normalized physical space (0-1)
const INITIAL_NODES = [
    {x: 0.1, y: 0.1},  // Node 1
    {x: 0.9, y: 0.1},  // Node 2
    {x: 0.9, y: 0.901},  // Node 3
    {x: 0.1, y: 0.9}   // Node 4
];

// Reference element nodes (ξ,η coordinates)
const REFERENCE_NODES = [
    {x: -1, y: -1}, // Node 1
    {x: 1, y: -1},  // Node 2
    {x: 1, y: 1},   // Node 3
    {x: -1, y: 1}   // Node 4
];

const FORCE_ARROW = {
    BASE_LENGTH: 100,
    MIN_LENGTH: 20,
    MAX_LENGTH: 500,
    HEAD: {
        LENGTH: 15,
        ANGLE: Math.PI / 6
    },
    HANDLE: {
        RADIUS: 5,
        HIT_RADIUS: 15
    },
    SCALE_RATIO: 0.25
};


// State
const state = {
    nodes: [...INITIAL_NODES.map(node => ({...node}))],
    selectedNode: null,
    isDragging: false,
    selectedShapeFunction: 0,
    showingDerivative: false,
    derivativeType: '',
    currentGaussPoint: { xi: 0, eta: 0 },
    currentGaussPointIndex: 0,
    isDraggingGaussPoint: false,
    standardGaussPoints: standardGaussPoints,
    refNodes: REFERENCE_NODES,
    view: 'U1',
    force: {
        length: FORCE_ARROW.BASE_LENGTH,
        angle: 0,
        isSelected: false,
        dragStart: null,
        getDisplayLength() {
            return this.length * FORCE_ARROW.SCALE_RATIO;
        }
    },
    // Add a place to store latest results if needed, or calculate on demand
    latestResults: {
        displacements: null,
        strainStress: null,
        nodalValues: {} // To store extrapolated nodal values
    }
};

// Canvas references
const canvases = {
    fem: null,
    reference: null,
    result: null,
    contexts: {
        fem: null,
        reference: null,
        result: null
    }
};

const controls = {
    config: {
        shapeFunctions: [
            { id: '1', label: 'N₁' , default: true},
            { id: '2', label: 'N₂' },
            { id: '3', label: 'N₃' },
            { id: '4', label: 'N₄' }
        ],
        derivatives: [
            { type: 'xi', label: '∂/∂ξ' },
            { type: 'eta', label: '∂/∂η' }
        ],
        views: [
            // Displacements
            { type: 'U1', label: 'U₁', group: 'Displacement', unit: 'µm', factor: 1e3, dot: 2, default: true},
            { type: 'U2', label: 'U₂', group: 'Displacement', unit: 'µm', factor: 1e3, dot: 2 },
            { type: 'Umag', label: '|U|', group: 'Displacement', unit: 'µm', factor: 1e3, dot: 2},
            // Strains
            { type: 'Exx', label: 'εxx', group: 'Strain', unit: 'µε', factor: 1e6, dot: 0 },
            { type: 'Eyy', label: 'εyy', group: 'Strain', unit: 'µε', factor: 1e6, dot: 0 },
            { type: 'Exy', label: 'γxy', group: 'Strain', unit: 'µrad', factor: 1e6, dot: 0 },
            // Stresses
            { type: 'Sxx', label: 'σxx', group: 'Stress', unit: 'MPa', factor: 1e-0, dot: 0 },
            { type: 'Syy', label: 'σyy', group: 'Stress', unit: 'MPa', factor: 1e-0, dot: 0 },
            { type: 'Sxy', label: 'τxy', group: 'Stress', unit: 'MPa', factor: 1e-0, dot: 0 },
        ]
    },

    handlers: {
        shapeFunction(index) {
            state.selectedShapeFunction = index;
            controls.updateButtons();
            renderer.reference();
            controls.updateReferenceTitle();
        },

        toggleDerivative() {
            state.showingDerivative = !state.showingDerivative;
            const derivativeControls = document.getElementById('derivativeControls');
            const toggleButton = document.querySelector('.toggle-button');

            toggleButton.textContent = state.showingDerivative ? 'Gradient' : 'Shape Function';
            toggleButton.classList.toggle('active', state.showingDerivative);
            derivativeControls.style.display = state.showingDerivative ? 'flex' : 'none';

            if (!state.showingDerivative) {
                state.derivativeType = '';
            } else {
                // Default to xi if enabling derivative view
                state.derivativeType = 'xi';
                document.querySelector('[data-type="xi"]').classList.add('selected');
            }

            renderer.reference();
            controls.updateReferenceTitle();
            controls.updateButtons();
        },

        derivative(type) {
            state.derivativeType = type;
            controls.updateButtons();
            renderer.reference();
            controls.updateReferenceTitle();
        },

        gaussPoint() {
            state.currentGaussPointIndex = (state.currentGaussPointIndex + 1) % state.standardGaussPoints.length;
            state.currentGaussPoint = {...state.standardGaussPoints[state.currentGaussPointIndex]};
            renderer.reference();
            updateLessons(state.currentGaussPoint, state.nodes, state.force);
        },

        viewChange(type) {
            state.view = type;
            controls.updateButtons(); // Use the new general updateButtons
            controls.updateResultTitle();
            // Re-render results with the new view
            // Calculation will happen within renderer.results
            renderer.results(); // Re-render results
        },

        // Add a helper to get current view config
        getCurrentViewConfig() {
            return controls.config.views.find(v => v.type === state.view);
        }
    },

    createButton(className, data, label, selected = false) {
        const button = document.createElement('button');
        button.className = className;
        button.textContent = label;
        if (data) {
            Object.entries(data).forEach(([key, value]) => {
                button.dataset[key] = value;
            });
        }
        if (selected) button.classList.add('selected');
        return button;
    },

    createButtonGroup(buttons, containerClass) {
        const group = document.createElement('div');
        group.className = containerClass;
        buttons.forEach(button => group.appendChild(button));
        return group;
    },

    createSeparator() {
        const separator = document.createElement('div');
        separator.className = 'control-separator';
        return separator;
    },

    updateButtons() {
        // Update shape function buttons
        document.querySelectorAll('.shape-button').forEach((button, index) => {
            button.classList.toggle('selected', index === state.selectedShapeFunction);
        });

        // Update derivative toggle and buttons
        const toggleButton = document.querySelector('.toggle-button');
         if (toggleButton) {
            toggleButton.textContent = state.showingDerivative ? 'Gradient' : 'Shape Function';
            toggleButton.classList.toggle('active', state.showingDerivative);
         }
        const derivativeControls = document.getElementById('derivativeControls');
         if (derivativeControls) {
            derivativeControls.style.display = state.showingDerivative ? 'flex' : 'none';
         }
        document.querySelectorAll('.derivative-button').forEach(button => {
            button.classList.toggle('selected', state.showingDerivative && button.dataset.type === state.derivativeType);
        });

        // Update view buttons (replaces displacement buttons update)
        document.querySelectorAll('.view-button').forEach(button => {
            button.classList.toggle('selected', button.dataset.type === state.view);
        });
    },

    updateReferenceTitle() {
        const title = document.getElementById('refElementTitle');
        if (!title) return;

        let titleText = state.showingDerivative
            ? `Shape function gradient ∂N${state.selectedShapeFunction + 1}/∂${state.derivativeType === 'xi' ? 'ξ' : 'η'}`
            : `Shape Function N${state.selectedShapeFunction + 1}`;

        title.innerHTML = `Reference Element (ξ,η)<br>${titleText}`;
    },

     updateResultTitle() {
        const title = document.getElementById('resultTitle');
        if (!title) return;
        const viewConfig = this.handlers.getCurrentViewConfig();
        if (viewConfig) {
            title.innerHTML = `Result<br>${viewConfig.group}: ${viewConfig.label} (${viewConfig.unit})`;
        } else {
             title.innerHTML = `Result<br>Unknown View`;
        }
    },

    initialize() {
        const buttonContainer = document.getElementById('controls');
        if (!buttonContainer) {
            console.error('Controls container not found');
            return;
        }
        buttonContainer.innerHTML = ''; // Clear existing

        // --- Row 1: Shape Functions & Derivatives ---
        const row1 = document.createElement('div');
        row1.className = 'control-row';

        // Shape Functions
        const shapeFunctionButtons = this.config.shapeFunctions.map(({ id, label, default: isDefault }) =>
            this.createButton('shape-button', { n: id }, label, isDefault && !state.showingDerivative) // Select default only if not showing derivatives initially
        );
        row1.appendChild(
            this.createButtonGroup(shapeFunctionButtons, 'button-group shape-functions')
        );
        row1.appendChild(this.createSeparator());

        // Toggle Button
        const toggleButton = this.createButton('toggle-button', null, 'Shape Function');
        row1.appendChild(toggleButton);
        row1.appendChild(this.createSeparator());

        // Derivative Controls (initially hidden)
        const derivativeContainer = document.createElement('div');
        derivativeContainer.id = 'derivativeControls';
        derivativeContainer.style.display = 'none';
        const derivativeButtons = this.config.derivatives.map(({ type, label }) =>
            this.createButton('derivative-button', { type }, label)
        );
        derivativeButtons.forEach(button => derivativeContainer.appendChild(button));
        row1.appendChild(derivativeContainer);
        row1.appendChild(this.createSeparator());

        // Gauss Button
        const gaussButton = this.createButton('gauss-button', null, 'Cycle Gauss Pt');
        row1.appendChild(gaussButton);

        buttonContainer.appendChild(row1);

        // --- Row 2: Result Views ---
        const row2 = document.createElement('div');
        row2.className = 'control-row';

        let currentGroup = '';
        this.config.views.forEach(({ type, label, group, default: isDefault }) => {
            if (group !== currentGroup) {
                if (currentGroup !== '') { // Add separator between groups
                   row2.appendChild(this.createSeparator());
                }
                // Optional: Add group label?
                const groupLabel = document.createElement('span');
                groupLabel.textContent = group + ': ';
                groupLabel.className = 'group-label';
                row2.appendChild(groupLabel);
                currentGroup = group;
            }
            // Add 'view-button' class for easier selection
            const button = this.createButton('view-button', { type }, label, isDefault);
            row2.appendChild(button);
        });
        buttonContainer.appendChild(row2);


        // Call initial updates
        this.updateButtons(); // Update selection states
        this.updateReferenceTitle();
        this.updateResultTitle();

        // Attach ALL event listeners
        this.attachEventListeners();
    },

    attachEventListeners() {
        // Shape function buttons
        document.querySelectorAll('.shape-button').forEach((button, index) => {
            button.addEventListener('click', () => this.handlers.shapeFunction(index));
        });

        // Toggle button
        const toggleButton = document.querySelector('.toggle-button');
        if (toggleButton) {
            toggleButton.addEventListener('click', () => this.handlers.toggleDerivative());
        }

        // Derivative buttons
        document.querySelectorAll('.derivative-button').forEach(button => {
            button.addEventListener('click', () =>
                this.handlers.derivative(button.dataset.type)
            );
        });

        // Gauss button
        const gaussButton = document.querySelector('.gauss-button');
        if (gaussButton) {
            gaussButton.addEventListener('click', () => this.handlers.gaussPoint());
        }

        // View buttons (replaces displacement buttons)
        document.querySelectorAll('.view-button').forEach(button => {
            button.addEventListener('click', () =>
                this.handlers.viewChange(button.dataset.type) // Use the renamed handler
            );
        });
    }
};

// Coordinate conversion functions
const coordConversion = {
    physicalToCanvas: (physicalPos, canvas) => ({
        x: physicalPos.x * canvas.width,
        y: physicalPos.y * canvas.height
    }),

    canvasToPhysical: (canvasPos, canvas) => ({
        x: canvasPos.x / canvas.width,
        y: canvasPos.y / canvas.height
    }),

    mapToReference: (x, y, canvas) => {
        const effectiveWidth = canvas.width - 2 * MARGIN;
        const effectiveHeight = canvas.height - 2 * MARGIN;
        return {
            xi: ((x - MARGIN) / effectiveWidth) * 2 - 1,
            eta: ((y - MARGIN) / effectiveHeight) * 2 - 1
        };
    },

    referenceToCanvas: (xi, eta, canvas) => {
        const effectiveWidth = canvas.width - 2 * MARGIN;
        const effectiveHeight = canvas.height - 2 * MARGIN;
        return {
            x: MARGIN + ((xi + 1) / 2) * effectiveWidth,
            y: MARGIN + ((eta + 1) / 2) * effectiveHeight
        };
    },

    findXiEta: (x, y, nodes, maxIterations = 10, tolerance = 1e-6) => {
        let xi = 0, eta = 0;

        for (let iter = 0; iter < maxIterations; iter++) {
            const N = calculateShapeFunctions(xi, eta);
            const dN = calculateDNdXi(xi, eta);

            let xCurrent = 0, yCurrent = 0;
            let dxdXi = 0, dxdEta = 0, dydXi = 0, dydEta = 0;

            for (let i = 0; i < 4; i++) {
                xCurrent += N[i] * nodes[i].x;
                yCurrent += N[i] * nodes[i].y;

                dxdXi += dN[i][0] * nodes[i].x;
                dxdEta += dN[i][1] * nodes[i].x;
                dydXi += dN[i][0] * nodes[i].y;
                dydEta += dN[i][1] * nodes[i].y;
            }

            const rx = x - xCurrent;
            const ry = y - yCurrent;
            const det = dxdXi * dydEta - dxdEta * dydXi;

            if (Math.abs(det) < 1e-10) return null;

            const dXi = (rx * dydEta - ry * dxdEta) / det;
            const dEta = (-rx * dydXi + ry * dxdXi) / det;

            xi += dXi;
            eta += dEta;

            if (Math.abs(dXi) < tolerance && Math.abs(dEta) < tolerance) {
                if (Math.abs(xi) <= 1.001 && Math.abs(eta) <= 1.001) {
                    return { xi, eta };
                }
                return null;
            }
        }
        return null;
    }
};

// Drawing functions
const renderer = {
    element() {
        const ctx = canvases.contexts.fem;
        const canvas = canvases.fem;
        if (!ctx || !canvas) return;

        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Draw element
        ctx.beginPath();
        const firstCanvasPos = coordConversion.physicalToCanvas(state.nodes[0], canvas);
        ctx.moveTo(firstCanvasPos.x, firstCanvasPos.y);

        for (let i = 1; i < state.nodes.length; i++) {
            const canvasPos = coordConversion.physicalToCanvas(state.nodes[i], canvas);
            ctx.lineTo(canvasPos.x, canvasPos.y);
        }
        ctx.closePath();
        ctx.strokeStyle = COLORS.element;
        ctx.stroke();

        // Draw nodes and their BCs
        state.nodes.forEach((node, i) => {
            const canvasPos = coordConversion.physicalToCanvas(node, canvas);

            // Draw node (similar to drawReferenceNodes)
            ctx.beginPath();
            ctx.arc(canvasPos.x, canvasPos.y, NODE_RADIUS, 0, Math.PI * 2);
            ctx.fillStyle = COLORS.node; // Keep original fill color, or change if needed
            ctx.fill();
            ctx.lineWidth = NODE_OUTLINE_WIDTH;
            ctx.strokeStyle = NODE_OUTLINE_COLOR;
            ctx.stroke();

            // Draw node label (similar to drawNodeLabel)
            ctx.font = NODE_LABEL_FONT;
            ctx.fillStyle = NODE_LABEL_COLOR; // Use specific color

            const text = `Node ${i + 1}`;
            const textMetrics = ctx.measureText(text);
            const textWidth = textMetrics.width;
            const textOffsetX = -textWidth / 2; // Center horizontally

            // Vertical offset logic (similar to reference, assuming nodes 0,1 are 'bottom' and 2,3 are 'top' in a quad)
            const textOffsetY = -(NODE_RADIUS + 5); // Above bottom nodes, below top nodes

            ctx.fillText(text, canvasPos.x + textOffsetX, canvasPos.y + textOffsetY);

            // Reset line width for subsequent drawings (like BCs, force)
            ctx.lineWidth = 1;

            // Draw fixed BC for node 4 (index 3)
            if (i === 3) {
                this.drawFixedBC(canvasPos);
            }

            // Draw rolling BC for node 3
            if (i === 2) {
                this.drawRollingBC(canvasPos);
            }
        });

        // Draw force arrow at node 2
        const node2Pos = coordConversion.physicalToCanvas(state.nodes[1], canvas);
        this.drawForceArrow(node2Pos, state.force);

    },

    drawFixedBC(pos) {
        const ctx = canvases.contexts.fem;
        const triangleSize = 15;

        ctx.beginPath();
        ctx.moveTo(pos.x - triangleSize, pos.y + 1.5 * triangleSize);
        ctx.lineTo(pos.x, pos.y); // Point at the node
        ctx.lineTo(pos.x + triangleSize, pos.y + 1.5 * triangleSize);
        ctx.closePath();
        ctx.strokeStyle = COLORS.boundaryCondition;
        ctx.stroke();
        ctx.fillStyle = COLORS.boundaryConditionFill;
        ctx.fill();
    },

    drawRollingBC(pos) {
        const ctx = canvases.contexts.fem;
        const circleRadius = 10;
        const offset = 15;

        ctx.beginPath();
        ctx.arc(pos.x, pos.y + offset, circleRadius, 0, Math.PI * 2);
        ctx.strokeStyle = COLORS.boundaryCondition;
        ctx.stroke();
        ctx.fillStyle = COLORS.boundaryConditionFill;
        ctx.fill();
    },

    drawForceArrow(pos, force) {
        const ctx = canvases.contexts.fem;
        const displayLength = force.getDisplayLength();

        // Calculate end point
        const endX = pos.x + displayLength * Math.cos(force.angle);
        const endY = pos.y + displayLength * Math.sin(force.angle);

        // Draw main line
        ctx.beginPath();
        ctx.moveTo(pos.x, pos.y);
        ctx.lineTo(endX, endY);
        ctx.strokeStyle = force.isSelected ? COLORS.forceSelected : COLORS.force;
        ctx.lineWidth = 2;
        ctx.stroke();

        // Draw arrowhead
        const headLength = FORCE_ARROW.HEAD.LENGTH;
        const headAngle = FORCE_ARROW.HEAD.ANGLE;

        ctx.beginPath();
        ctx.moveTo(endX, endY);
        ctx.lineTo(
            endX - headLength * Math.cos(force.angle - headAngle),
            endY - headLength * Math.sin(force.angle - headAngle)
        );
        ctx.moveTo(endX, endY);
        ctx.lineTo(
            endX - headLength * Math.cos(force.angle + headAngle),
            endY - headLength * Math.sin(force.angle + headAngle)
        );
        ctx.stroke();

        // Draw handle
        ctx.beginPath();
        ctx.arc(endX, endY, FORCE_ARROW.HANDLE.RADIUS, 0, Math.PI * 2);
        ctx.fillStyle = force.isSelected ? COLORS.forceSelected : COLORS.force;
        ctx.fill();

        // Draw magnitude
        ctx.fillStyle = COLORS.text;
        ctx.font = '12px Arial';
        const magnitude = Math.round(force.length);
        ctx.fillText(`${magnitude} N`, endX + 10, endY + 10);

        ctx.lineWidth = 1;
    },

    reference() {
        const ctx = canvases.contexts.reference;
        const canvas = canvases.reference;
        if (!ctx || !canvas) return;

        ctx.clearRect(0, 0, canvas.width, canvas.height);

        const effectiveWidth = canvas.width - 2 * MARGIN;
        const effectiveHeight = canvas.height - 2 * MARGIN;

        // First pass: find min and max values
        let minValue = Infinity;
        let maxValue = -Infinity;
        const values = new Array(effectiveWidth * effectiveHeight);

        for(let y = 0; y < effectiveHeight; y++) {
            for(let x = 0; x < effectiveWidth; x++) {
                const xi = (x / effectiveWidth) * 2 - 1;
                const eta = (y / effectiveHeight) * 2 - 1;

                let value;
                if (state.showingDerivative) {
                    const derivatives = calculateDNdXi(xi, eta);
                    value = state.derivativeType === 'xi' ?
                        derivatives[state.selectedShapeFunction][0] :
                        derivatives[state.selectedShapeFunction][1];
                } else {
                    const N = calculateShapeFunctions(xi, eta);
                    value = N[state.selectedShapeFunction];
                }

                const index = y * effectiveWidth + x;
                values[index] = value;
                minValue = Math.min(minValue, value);
                maxValue = Math.max(maxValue, value);
            }
        }

        // Create and draw image data
        const pixelData = new Uint8ClampedArray(effectiveWidth * effectiveHeight * 4);
        const valueRange = maxValue - minValue;

        for(let y = 0; y < effectiveHeight; y++) {
            for(let x = 0; x < effectiveWidth; x++) {
                const index = y * effectiveWidth + x;
                const value = values[index];
                const normalizedValue = valueRange === 0 ? 0.5 : (value - minValue) / valueRange;
                const color = this.getDisplacementColor(normalizedValue);

                const pixelIndex = index * 4;
                pixelData[pixelIndex] = color.r;
                pixelData[pixelIndex + 1] = color.g;
                pixelData[pixelIndex + 2] = color.b;
                pixelData[pixelIndex + 3] = 255;
            }
        }

        const imageData = new ImageData(pixelData, effectiveWidth, effectiveHeight);
        ctx.putImageData(imageData, MARGIN, MARGIN);

        // Draw axes
        this.drawReferenceAxes(ctx, canvas);

        // Draw Gauss points
        this.drawGaussPoints(ctx, canvas);

        // Draw reference element outline
        this.drawReferenceOutline(ctx, canvas);

        // Draw nodes
        this.drawReferenceNodes(ctx, canvas);

        // Draw Reference Legend (similar to displacement legend)
        this.drawLegend(ctx, canvas, minValue, maxValue, 1);
    },

    drawReferenceAxes(ctx, canvas) {
        ctx.beginPath();
        ctx.moveTo(MARGIN, canvas.height/2);
        ctx.lineTo(canvas.width - MARGIN, canvas.height/2);
        ctx.moveTo(canvas.width/2, MARGIN);
        ctx.lineTo(canvas.width/2, canvas.height - MARGIN);
        ctx.strokeStyle = 'rgba(128,128,128,0.5)';
        ctx.stroke();
    },

    drawGaussPoints(ctx, canvas) {
        // Draw standard Gauss points
        state.standardGaussPoints.forEach(point => {
            const canvasPoint = coordConversion.referenceToCanvas(point.xi, point.eta, canvas);
            ctx.beginPath();
            ctx.arc(canvasPoint.x, canvasPoint.y, 4, 0, Math.PI * 2);
            ctx.fillStyle = 'rgba(255, 255, 255, 1.0)';
            ctx.fill();
            ctx.strokeStyle = 'black';
            ctx.stroke();
        });

        // Draw current Gauss point
        const gaussCanvasPoint = coordConversion.referenceToCanvas(
            state.currentGaussPoint.xi,
            state.currentGaussPoint.eta,
            canvas
        );

        // Draw crosshair and point
        this.drawGaussPointCrosshair(ctx, gaussCanvasPoint);
    },

    drawGaussPointCrosshair(ctx, point) {
        ctx.beginPath();
        ctx.arc(point.x, point.y, 4, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(255, 0, 255, 0.5)';
        ctx.fill();
        ctx.strokeStyle = 'purple';
        ctx.lineWidth = 2;
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(point.x - 10, point.y);
        ctx.lineTo(point.x + 10, point.y);
        ctx.moveTo(point.x, point.y - 10);
        ctx.lineTo(point.x, point.y + 10);
        ctx.stroke();
        ctx.lineWidth = 1;
    },

    drawReferenceOutline(ctx, canvas) {
        ctx.beginPath();
        const startPoint = coordConversion.referenceToCanvas(
            state.refNodes[0].x,
            state.refNodes[0].y,
            canvas
        );
        ctx.moveTo(startPoint.x, startPoint.y);

        for(let i = 1; i < state.refNodes.length; i++) {
            const point = coordConversion.referenceToCanvas(
                state.refNodes[i].x,
                state.refNodes[i].y,
                canvas
            );
            ctx.lineTo(point.x, point.y);
        }
        ctx.closePath();
        ctx.strokeStyle = 'blue';
        ctx.stroke();
    },

    drawReferenceNodes(ctx, canvas) {
        state.refNodes.forEach((node, i) => {
            const point = coordConversion.referenceToCanvas(node.x, node.y, canvas);
            const nodeRadius = 8;

            ctx.beginPath();
            ctx.arc(point.x, point.y, nodeRadius, 0, Math.PI * 2);
            ctx.fillStyle = i === state.selectedShapeFunction ? 'yellow' : 'red';
            ctx.fill();
            ctx.lineWidth = 2;
            ctx.strokeStyle = 'black';
            ctx.stroke();

            this.drawNodeLabel(ctx, point, i, nodeRadius);
            ctx.lineWidth = 1;
        });
    },

    drawNodeLabel(ctx, point, index, nodeRadius) {
        ctx.font = '14px Arial';
        ctx.fillStyle = 'black';

        const text = `Node ${index + 1}`;
        const textWidth = ctx.measureText(text).width;
        const textOffsetX = -textWidth / 2;
        const textOffsetY = (index === 0 || index === 1) ?
            -(nodeRadius + 5) :
            nodeRadius + 15;

        ctx.fillText(text, point.x + textOffsetX, point.y + textOffsetY);
    },

     drawLegend(ctx, canvas, minValue, maxValue, Dot) {
        // --- Box Configuration ---
        const boxPadding = 5; // Padding inside and outside the box
        const numDigits = Math.max(
            Math.abs(minValue) > 1e-3 ? Math.max(1, Math.floor(Math.log10(Math.abs(minValue)))) + 1 : 0,
            Math.abs(maxValue) > 1e-3 ? Math.max(1, Math.floor(Math.log10(Math.abs(maxValue)))) + 1 : 0,
            0 // ensure at least 0
        );
        const estimatedMaxTextWidth = (numDigits + Dot + (Dot > 0 ? 1 : 0) + 2) * 6 + 10; // Rough estimate (adjust multiplier as needed)
        const rectWidth = 10;
        const rectHeight = 5;
        const textPadding = 5; // Space between text and color rectangle
        const lineSpacing = 5; // Space between max and min rows

        // --- Determine Content Dimensions ---
        const contentWidth = estimatedMaxTextWidth + textPadding + rectWidth;
        let contentHeight;
        const isSingleValue = Math.abs(maxValue - minValue) <= 1e-9 * (Math.abs(minValue) + Math.abs(maxValue) + 1) || !isFinite(minValue) || !isFinite(maxValue);


        if (isSingleValue) {
            if (!isFinite(minValue)) {
                console.warn("Legend: Cannot draw legend for non-finite values.");
                return; // Don't draw anything if the single value isn't finite
            }
            contentHeight = rectHeight; // Only one row needed
        } else {
            contentHeight = rectHeight + lineSpacing + rectHeight; // Height for two rows (max and min)
        }

        // --- Calculate Box Dimensions & Position ---
        const boxWidth = contentWidth + 2 * boxPadding;
        const boxHeight = contentHeight + 2 * boxPadding;
        const boxX = boxPadding; // Position box near left edge
        const boxY = canvas.height / 2 - boxHeight / 2; // Center box vertically

        // --- Draw the Background Box ---
        ctx.save(); // Save context state before drawing the box
        ctx.fillStyle = 'white';
        ctx.fillRect(boxX, boxY, boxWidth, boxHeight);
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1;
        ctx.strokeRect(boxX, boxY, boxWidth, boxHeight);
        ctx.restore(); // Restore context state

        // --- Calculate Content Area Origin ---
        const contentOriginX = boxX + boxPadding;
        const contentOriginY = boxY + boxPadding;
        const contentRightEdge = boxX + boxWidth - boxPadding;

        // --- Prepare for Drawing Legend Items ---
        ctx.font = '10px Arial';
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle'; // Align text vertically with rectangle center

        // --- Draw Legend Content ---
        const rectX = contentRightEdge - rectWidth;
        const textX = contentRightEdge - rectWidth - textPadding;

        const formatNumber = (num, precision) => {
             // Handle very small or very large numbers with exponential notation
             //if (Math.abs(num) > 1e5 || (Math.abs(num) < 1e-3 && num !== 0)) {
             //    return num.toExponential(precision > 0 ? precision -1 : 0);
             //}
             return num.toFixed(precision);
         };


        if (isSingleValue) {
            const midColor = this.getDisplacementColor(0.5); // Or 0 or 1
            const midText = formatNumber(minValue, Dot); // Use formatter and Dot
            const itemY = contentOriginY + contentHeight / 2 - rectHeight / 2; // Center vertically
            const textY = itemY + rectHeight / 2;

            ctx.fillStyle = 'black';
            ctx.fillText(midText, textX, textY);
            ctx.fillStyle = `rgb(${midColor.r}, ${midColor.g}, ${midColor.b})`;
            ctx.fillRect(rectX, itemY, rectWidth, rectHeight);
            ctx.strokeStyle = 'black';
            ctx.strokeRect(rectX, itemY, rectWidth, rectHeight);

        } else {
            const minColor = this.getDisplacementColor(0);
            const maxColor = this.getDisplacementColor(1);
            const minText = formatNumber(minValue, Dot); // Use formatter and Dot
            const maxText = formatNumber(maxValue, Dot); // Use formatter and Dot

            // Max Value (Top Row)
            const maxItemY = contentOriginY;
            const maxTextY = maxItemY + rectHeight / 2;
            ctx.fillStyle = 'black';
            ctx.fillText(maxText, textX, maxTextY);
            ctx.fillStyle = `rgb(${maxColor.r}, ${maxColor.g}, ${maxColor.b})`;
            ctx.fillRect(rectX, maxItemY, rectWidth, rectHeight);
            ctx.strokeStyle = 'black';
            ctx.strokeRect(rectX, maxItemY, rectWidth, rectHeight);

            // Min Value (Bottom Row)
            const minItemY = maxItemY + rectHeight + lineSpacing;
            const minTextY = minItemY + rectHeight / 2;
            ctx.fillStyle = 'black';
            ctx.fillText(minText, textX, minTextY);
            ctx.fillStyle = `rgb(${minColor.r}, ${minColor.g}, ${minColor.b})`;
            ctx.fillRect(rectX, minItemY, rectWidth, rectHeight);
            ctx.strokeStyle = 'black';
            ctx.strokeRect(rectX, minItemY, rectWidth, rectHeight);
        }

        // Reset context properties
        ctx.textAlign = 'left';
        ctx.textBaseline = 'alphabetic';
        ctx.lineWidth = 1;
    },

    results() {
        const ctx = canvases.contexts.result;
        const canvas = canvases.result;
        if (!ctx || !canvas) return;

        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // --- 1. Calculate Displacements and Strain/Stress ---
        const displacements = calculateDisplacements(state.nodes, state.force);
        // Ensure displacements is always a valid array, even if calculation fails
         if (!displacements || displacements.length !== 4) {
             console.error("Failed to calculate displacements.");
             // Optionally draw meshes without contours or return
              this.drawUndeformedMesh();
              this.drawDeformedMesh(state.nodes); // Draw undeformed as deformed
             return;
         }
        const strainStressResults = calculateStrainAndStress(state.nodes, displacements);
         if (!strainStressResults || strainStressResults.length !== 4) {
             console.error("Failed to calculate strain/stress.");
             // Optionally draw meshes without contours or return
              this.drawUndeformedMesh();
              this.drawDeformedMesh(state.nodes); // Draw undeformed as deformed
              return;
         }

        // Store for potential later use (e.g., Gauss point display)
        state.latestResults.displacements = displacements;
        state.latestResults.strainStress = strainStressResults;


        // --- 2. Determine Nodal Values for Current View ---
        const currentView = state.view;
        let nodalValues;
        let componentKey = null; // For extrapolation

        if (currentView.startsWith('U')) {
            if (currentView === 'U1') nodalValues = displacements.map(d => d.x);
            else if (currentView === 'U2') nodalValues = displacements.map(d => d.y);
            else { // Umag
                nodalValues = displacements.map(d => Math.sqrt(d.x*d.x + d.y*d.y));
            }
        } else {
            // Extrapolate Stress/Strain
             componentKey = currentView; // e.g., 'Sxx', 'Eyy'
             nodalValues = extrapolateGaussToNodes(strainStressResults, componentKey);
             state.latestResults.nodalValues[componentKey] = nodalValues; // Cache if needed
        }

        // --- 3. Prepare Deformed Geometry ---
        const scaledNodes = state.nodes.map((node, i) => ({
            x: node.x + displacements[i].x * DISPLACEMENT_SCALE_FACTOR,
            y: node.y + displacements[i].y * DISPLACEMENT_SCALE_FACTOR
        }));
        const canvasNodes = scaledNodes.map(node =>
            coordConversion.physicalToCanvas(node, canvas)
        );

        // --- 4. Render Contour Plot (Pixel Loop) ---
        const minX = Math.min(...canvasNodes.map(n => n.x));
        const maxX = Math.max(...canvasNodes.map(n => n.x));
        const minY = Math.min(...canvasNodes.map(n => n.y));
        const maxY = Math.max(...canvasNodes.map(n => n.y));

        const width = Math.ceil(maxX - minX);
        const height = Math.ceil(maxY - minY);

        // Check for zero width/height (can happen if element is degenerate)
         if (width <= 0 || height <= 0) {
             console.warn("Degenerate element geometry, skipping contour plot.");
              this.drawUndeformedMesh();
              this.drawDeformedMesh(scaledNodes); // Still draw deformed shape
             return;
         }


        const pixelData = new Uint8ClampedArray(width * height * 4);
        let minValue = Infinity;
        let maxValue = -Infinity;
        const values = new Array(width * height).fill(undefined); // Use undefined marker

        // Find Min/Max by interpolating at nodes first (faster approximation)
         // or use the actual nodalValues directly
         minValue = Math.min(...nodalValues);
         maxValue = Math.max(...nodalValues);
         // Refine min/max in the loop if needed, but nodal values are often sufficient

        for(let y = 0; y < height; y++) {
            for(let x = 0; x < width; x++) {
                const canvasX = x + minX;
                const canvasY = y + minY;

                // Check if pixel is roughly inside the element polygon (optimization)
                // Optional, but can speed up large canvases
                // if (!isPointInPolygon({x: canvasX, y: canvasY}, canvasNodes)) continue;

                const physX = canvasX / canvas.width;
                const physY = canvasY / canvas.height;

                // Map canvas pixel to xi, eta in DEFORMED element
                const result = coordConversion.findXiEta(physX, physY, scaledNodes);

                if (result) {
                    const { xi, eta } = result;
                    // Ensure xi, eta are within bounds (or close enough) after inversion
                    if (Math.abs(xi) <= 1.01 && Math.abs(eta) <= 1.01) {
                        const N = calculateShapeFunctions(xi, eta);
                        let value = 0;
                        for(let i = 0; i < 4; i++) {
                            value += N[i] * nodalValues[i];
                        }

                        const index = y * width + x;
                        values[index] = value; // Store value
                        // Update min/max if needed (safer than just using nodal extremes)
                        minValue = Math.min(minValue, value);
                        maxValue = Math.max(maxValue, value);
                    }
                }
            }
        }

        // Second pass: create color data
        const valueRange = maxValue - minValue;
        for(let y = 0; y < height; y++) {
            for(let x = 0; x < width; x++) {
                const index = y * width + x;
                const value = values[index]; // Retrieve stored value

                if (value !== undefined) { // Check if value was calculated
                    const normalizedValue = valueRange === 0 ? 0.5 : (value - minValue) / valueRange;
                    const color = this.getDisplacementColor(normalizedValue); // Reusing color scheme

                    const pixelIndex = index * 4;
                    pixelData[pixelIndex] = color.r;
                    pixelData[pixelIndex + 1] = color.g;
                    pixelData[pixelIndex + 2] = color.b;
                    pixelData[pixelIndex + 3] = 255; // Fully opaque
                }
                 // else: leave pixel transparent (or set a background color)
            }
        }

        // Create and draw the image data
        const imageData = new ImageData(pixelData, width, height);
        ctx.putImageData(imageData, minX, minY);

        // --- 5. Draw Overlays ---
        this.drawUndeformedMesh();
        this.drawDeformedMesh(scaledNodes); // Use the calculated scaledNodes
        this.drawDisplacementScale(); // Keep this relevant

        // --- 6. Draw Legend ---
        const viewConfig = controls.handlers.getCurrentViewConfig();
        if (viewConfig) {
            // Apply scaling factor for display units (e.g., Pa -> MPa)
            this.drawLegend(ctx, canvas, minValue * viewConfig.factor, maxValue * viewConfig.factor, viewConfig.dot);
        }

        // --- 7. Draw Gauss Point Values (Optional Enhancement) ---
         this.drawGaussPointValues(ctx, canvas, scaledNodes, state.latestResults.strainStress, viewConfig);


    }, // End of renderer.results

     // New function to draw values at Gauss Points on the result canvas
    drawGaussPointValues(ctx, canvas, scaledNodes, strainStressResults, viewConfig) {
        // Basic validation for inputs
        if (!ctx || !canvas || !scaledNodes || scaledNodes.length < 4 || !strainStressResults || !viewConfig) {
            // console.warn("drawGaussPointValues: Missing required arguments.");
            return;
        }

        ctx.font = '10px Arial';
        ctx.fillStyle = 'black';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';

        // Iterate through the results array directly
        strainStressResults.forEach((result, index) => {
            // --- Robustness Check ---
            if (!result) {
                // console.warn(`[GP ${index + 1}] Skipping draw: No result data found.`);
                return; // Skip this iteration if result is missing
            }

            // Use the gaussPoint data stored within the result if available, otherwise fallback
            // This ensures we use the correct coordinates even if it's an error result
            const gp = result.gaussPoint || state.standardGaussPoints[index];
            if (!gp || typeof gp.xi !== 'number' || typeof gp.eta !== 'number') {
                console.warn(`[GP ${index + 1}] Skipping draw: Invalid or missing Gauss point coordinates.`);
                return; // Cannot calculate position without valid GP coords
            }

            let canvasPos;
            try {
                // Calculate deformed physical coordinates of the Gauss Point
                const N = calculateShapeFunctions(gp.xi, gp.eta);
                let deformedPhysX = 0;
                let deformedPhysY = 0;
                // Ensure we don't go out of bounds for scaledNodes
                for (let i = 0; i < Math.min(4, scaledNodes.length); i++) {
                    if (scaledNodes[i]) { // Check node exists
                        deformedPhysX += N[i] * scaledNodes[i].x;
                        deformedPhysY += N[i] * scaledNodes[i].y;
                    } else {
                         console.warn(`[GP ${index + 1}] Skipping draw: Missing node data at index ${i}.`);
                         return; // Cannot calculate position accurately
                    }
                }
                // Convert deformed physical coordinates to canvas coordinates
                canvasPos = coordConversion.physicalToCanvas({ x: deformedPhysX, y: deformedPhysY }, canvas);

            } catch (e) {
                console.error(`[GP ${index + 1}] Error calculating position:`, e);
                return; // Skip if position calculation fails
            }

            // --- Check for Calculation Error in Result ---
            if (result.error) {
                // Optionally, draw an indicator that calculation failed for this point
                try {
                    ctx.save(); // Save context state
                    ctx.fillStyle = 'red';
                    ctx.font = 'bold 12px Arial';
                    ctx.fillText('!', canvasPos.x, canvasPos.y + 1); // Draw red exclamation mark
                    ctx.restore(); // Restore context state (font, fillStyle)
                    // console.warn(`[GP ${index + 1}] Calculation error reported: ${result.error}`);
                } catch(drawError) {
                     console.error(`[GP ${index + 1}] Error drawing error marker:`, drawError);
                }
                return; // Skip drawing the value due to the reported error
            }

            // --- If no error, proceed to get and draw the value ---
            let value;
            try {
                // Check if the required data structure exists before accessing it
                if (viewConfig.group === 'Stress' && result.stress) {
                     const key = viewConfig.type === 'Sxx' ? 'σxx' : viewConfig.type === 'Syy' ? 'σyy' : 'τxy';
                     value = result.stress[key];
                 } else if (viewConfig.group === 'Strain' && result.strain) {
                      const key = viewConfig.type === 'Exx' ? 'εxx' : viewConfig.type === 'Eyy' ? 'εyy' : 'γxy';
                     value = result.strain[key];
                 } else if (viewConfig.group === 'Displacement') {
                     return; // Don't draw values for displacement at GPs for now
                 } else {
                     // Handle cases where the group is Stress/Strain but the object is unexpectedly missing
                     // console.warn(`[GP ${index + 1}] Missing ${viewConfig.group} data in valid result.`);
                     return;
                 }

                 // Check if a valid numeric value was retrieved
                if (typeof value === 'number' && isFinite(value)) {
                    const displayValue = (value * viewConfig.factor).toFixed(viewConfig.dot);

                     // Draw a small background box for readability
                     const textWidth = ctx.measureText(displayValue).width;
                     ctx.fillStyle = 'rgba(255, 255, 255, 0.7)'; // Semi-transparent white
                     ctx.fillRect(canvasPos.x - textWidth / 2 - 2, canvasPos.y - 7, textWidth + 4, 14);

                     // Draw the text value
                     ctx.fillStyle = 'black'; // Ensure text color is reset
                    ctx.fillText(displayValue, canvasPos.x, canvasPos.y);
                 } else {
                     // Handle cases where the key exists but the value is not a valid number (e.g., NaN)
                     // console.warn(`[GP ${index + 1}] Invalid value found for ${viewConfig.group}/${viewConfig.type}:`, value);
                     // Optionally draw 'NaN' or similar
                     ctx.save();
                     ctx.fillStyle = 'rgba(255, 255, 255, 0.7)';
                     ctx.fillRect(canvasPos.x - 10, canvasPos.y - 7, 20, 14);
                     ctx.fillStyle = 'orange'; // Use a different color for NaN/invalid
                     ctx.fillText('N/A', canvasPos.x, canvasPos.y);
                     ctx.restore();
                 }

            } catch (e) {
                 // This catch now handles errors *after* the initial error check
                 // (e.g., issue with viewConfig, key access on potentially malformed but not 'error' result, drawing errors)
                 console.error(`[GP ${index + 1}] Error processing or drawing value:`, e, "Result object:", result);
                 // Optionally draw an error marker here too
                 try {
                     ctx.save();
                     ctx.fillStyle = 'red';
                     ctx.font = 'bold 12px Arial';
                     ctx.fillText('!', canvasPos.x, canvasPos.y + 1);
                     ctx.restore();
                 } catch(drawError) {
                     console.error(`[GP ${index + 1}] Error drawing error marker in catch block:`, drawError);
                 }
            }
        });

        // Reset alignment
         ctx.textAlign = 'left';
         ctx.textBaseline = 'alphabetic';
    },

    getDisplacementColor(normalizedValue) {
        const spectralColors = [
            {r: 94, g: 79, b: 162},   // #5e4fa2
            {r: 50, g: 136, b: 189},  // #3288bd
            {r: 102, g: 194, b: 165}, // #66c2a5
            {r: 171, g: 221, b: 164}, // #abdda4
            {r: 230, g: 245, b: 152}, // #e6f598
            {r: 255, g: 255, b: 191}, // #ffffbf
            {r: 254, g: 224, b: 139}, // #fee08b
            {r: 253, g: 174, b: 97},  // #fdae61
            {r: 244, g: 109, b: 67},  // #f46d43
            {r: 213, g: 62, b: 79},   // #d53e4f
            {r: 158, g: 1, b: 66}     // #9e0142
        ];

        // Ensure the normalized value is between 0 and 1
        const clampedValue = Math.max(0, Math.min(1, normalizedValue));

        // Simply select the color based on the normalized value
        const index = Math.min(
            Math.floor(clampedValue * spectralColors.length),
            spectralColors.length - 1
        );

        return spectralColors[index];
    },

    drawUndeformedMesh() {
        const ctx = canvases.contexts.result;
        const canvas = canvases.result;

        ctx.beginPath();
        const firstNode = coordConversion.physicalToCanvas(state.nodes[0], canvas);
        ctx.moveTo(firstNode.x, firstNode.y);

        for(let i = 1; i < state.nodes.length; i++) {
            const node = coordConversion.physicalToCanvas(state.nodes[i], canvas);
            ctx.lineTo(node.x, node.y);
        }
        ctx.closePath();
        ctx.setLineDash([5, 5]);
        ctx.strokeStyle = 'gray';
        ctx.stroke();
        ctx.setLineDash([]);
    },

    drawDeformedMesh(scaledNodes) {
        const ctx = canvases.contexts.result;
        const canvas = canvases.result;

        ctx.beginPath();
        const firstNode = coordConversion.physicalToCanvas(scaledNodes[0], canvas);
        ctx.moveTo(firstNode.x, firstNode.y);

        for(let i = 1; i < scaledNodes.length; i++) {
            const node = coordConversion.physicalToCanvas(scaledNodes[i], canvas);
            ctx.lineTo(node.x, node.y);
        }
        ctx.closePath();
        ctx.strokeStyle = 'blue';
        ctx.lineWidth = 2;
        ctx.stroke();
        ctx.lineWidth = 1;
    },

    drawDisplacementScale() {
        const ctx = canvases.contexts.result;
        ctx.fillStyle = 'black';
        ctx.font = '12px Arial';
        ctx.fillText(`Displacement Scale: ${DISPLACEMENT_SCALE_FACTOR}×`, 10, 20);
    }
};

// Event handlers
const handlers = {
    // Renamed from getMousePos to handle both mouse and touch
    getEventPos(canvas, evt) {
        const rect = canvas.getBoundingClientRect();
        let touch = null;

        if (evt.type.startsWith('touch')) {
            // For touchstart, touchend, touchcancel use the relevant touch from changedTouches
            // For touchmove, use the first active touch from touches
            touch = (evt.type === 'touchmove') ? evt.touches[0] : evt.changedTouches[0];
        }

        const clientX = touch ? touch.clientX : evt.clientX;
        const clientY = touch ? touch.clientY : evt.clientY;

        // Check if coordinates are valid before calculating relative position
        if (clientX === undefined || clientY === undefined) {
             console.warn("Could not determine event coordinates.");
             return { x: 0, y: 0 }; // Return origin as fallback
        }

        return {
            x: clientX - rect.left,
            y: clientY - rect.top
        };
    },

    GlobalMouseUp(event) {
        const wasDragging = state.isDragging || state.isDraggingGaussPoint;

        state.isDragging = false;
        state.isDraggingGaussPoint = false;
        state.selectedNode = null;
        state.force.isSelected = false;
        state.force.dragStart = null;

         if (wasDragging) {
            // console.log("Global up detected, drag ended.");
         }
    },

    FemMouseDown(event) {
        const pos = this.getEventPos(canvases.fem, event);
        const nodeIndex = this.findSelectedNode(pos, canvases.fem); // Find *any* node

        if (this.isForceArrowSelected(pos)) {
            this.initializeForceArrowDrag(pos);
        } else if (nodeIndex !== -1) { // Check if *any* node was found
            this.initializeNodeDrag(nodeIndex); // Allow dragging any found node
        }
        // No preventDefault here, allow default actions unless dragging starts
    },

    FemMouseMove(event) {
        // Only proceed if a drag operation is active (either node or force)
        if (!state.isDragging) return;

        // Prevent scrolling/zooming on touch devices during drag
        if (event.type === 'touchmove') {
            event.preventDefault();
        }

        const pos = this.getEventPos(canvases.fem, event);

        if (state.force.isSelected) {
            this.updateForceArrow(pos);
        } else if (state.selectedNode !== null) { // Check if a node is selected
            this.updateNodePosition(pos); // Update the selected node's position
        }

        // Only redraw/recalculate if something was actually dragged
        renderer.element();
        renderer.results();
        updateLessons(state.currentGaussPoint, state.nodes, state.force);
    },

     FemMouseUp(event) {
        // GlobalMouseUp handles the state reset, this handler might be useful
        // for snapping or other finalization specific to the FEM canvas,
        // but currently not needed for node/force dragging.
        if (state.isDragging) {
            // console.log("FEM Up/End");
        }
        // Reset flags are handled globally by GlobalMouseUp
    },

    RefMouseDown(event) {
        const pos = this.getEventPos(canvases.reference, event);
        const gaussCanvasPoint = coordConversion.referenceToCanvas(
            state.currentGaussPoint.xi,
            state.currentGaussPoint.eta,
            canvases.reference
        );

        // Use a slightly larger radius for touch targets
        const hitRadius = (event.type.startsWith('touch')) ? 20 : 15;
        const dx = pos.x - gaussCanvasPoint.x;
        const dy = pos.y - gaussCanvasPoint.y;

        if (dx * dx + dy * dy < hitRadius * hitRadius) {
            state.isDraggingGaussPoint = true;
            // No preventDefault here, wait for move to confirm drag intention
        }
    },

    RefMouseMove(event) {
        if (state.isDraggingGaussPoint) {
             // Prevent scrolling/zooming on touch devices during drag
             if (event.type === 'touchmove') {
                 event.preventDefault();
             }

            const pos = this.getEventPos(canvases.reference, event);

            // Clamp mouse position within effective canvas area
            const clampedX = Math.max(MARGIN, Math.min(pos.x, canvases.reference.width - MARGIN));
            const clampedY = Math.max(MARGIN, Math.min(pos.y, canvases.reference.height - MARGIN));

            // Convert to reference coordinates
            const referenceCoords = coordConversion.mapToReference(
                clampedX,
                clampedY,
                canvases.reference
            );

            // Update Gauss point position
            state.currentGaussPoint.xi = referenceCoords.xi;
            state.currentGaussPoint.eta = referenceCoords.eta;
             // Reset index when manually dragging away from a standard point
             state.currentGaussPointIndex = null; // Or -1, depending on convention

            // Redraw and update lessons
            renderer.reference();
            updateLessons(state.currentGaussPoint, state.nodes, state.force);
        }
    },

    RefMouseUp(event) {
        if (state.isDraggingGaussPoint) {
            // Reset drag state immediately
            state.isDraggingGaussPoint = false;
            // console.log("Ref Up/End");

            // --- Snap Logic ---
            const threshold = 0.2; // Snap distance in reference coordinates
            let closestPointInfo = null;
            let minDistance = Infinity;

            state.standardGaussPoints.forEach((point, index) => {
                const dx = point.xi - state.currentGaussPoint.xi;
                const dy = point.eta - state.currentGaussPoint.eta;
                const distance = Math.sqrt(dx * dx + dy * dy);

                if (distance < minDistance) {
                    minDistance = distance;
                    closestPointInfo = { point, index };
                }
            });

            // Snap to closest standard point if within threshold
            if (closestPointInfo && minDistance < threshold) {
                state.currentGaussPoint.xi = closestPointInfo.point.xi;
                state.currentGaussPoint.eta = closestPointInfo.point.eta;
                state.currentGaussPointIndex = closestPointInfo.index; // Store the index of the snapped point
                // Redraw and update after snapping
                renderer.reference();
                updateLessons(state.currentGaussPoint, state.nodes, state.force);
            }
            // --- End Snap Logic ---
        }
         // Reset flags are handled globally by GlobalMouseUp
    },

    // --- Helper Methods ---

    findSelectedNode(pos, canvas) {
        // Increase hit radius slightly for touch
        // NODE_RADIUS should be defined globally or passed in
        const effectiveNodeRadius = (event.type && event.type.startsWith('touch'))
            ? NODE_RADIUS * 1.5
            : NODE_RADIUS;

        return state.nodes.findIndex(node => { // No filter based on index needed here
            const canvasPos = coordConversion.physicalToCanvas(node, canvas);
            const distance = Math.hypot(canvasPos.x - pos.x, canvasPos.y - pos.y);
            return distance < effectiveNodeRadius;
        });
    },

    isForceArrowSelected(pos) {
        const node1Canvas = coordConversion.physicalToCanvas(state.nodes[1], canvases.fem);
        const arrowEndX = node1Canvas.x + state.force.getDisplayLength() * Math.cos(state.force.angle);
        const arrowEndY = node1Canvas.y + state.force.getDisplayLength() * Math.sin(state.force.angle);

        const distToArrowEnd = Math.hypot(pos.x - arrowEndX, pos.y - arrowEndY);
        // FORCE_ARROW.HANDLE.HIT_RADIUS should be defined globally or passed in
        return distToArrowEnd < FORCE_ARROW.HANDLE.HIT_RADIUS;
    },

    initializeForceArrowDrag(pos) {
        state.force.isSelected = true;
        state.force.dragStart = {
            x: pos.x,
            y: pos.y,
            angle: state.force.angle,
            length: state.force.length
        };
        state.isDragging = true; // Set the general dragging flag
        // console.log("Dragging force arrow started");
    },

    initializeNodeDrag(nodeIndex) {
        // No need to check nodeIndex against 2 or 3 here anymore
        state.selectedNode = nodeIndex;
        state.isDragging = true; // Set the general dragging flag
        // console.log(`Dragging node ${nodeIndex + 1} started`);
    },

    updateForceArrow(pos) {
        const node1Canvas = coordConversion.physicalToCanvas(state.nodes[1], canvases.fem);
        const dx = pos.x - node1Canvas.x;
        const dy = pos.y - node1Canvas.y;

        // Update angle
        state.force.angle = Math.atan2(dy, dx);

        // Calculate and constrain length
        // FORCE_ARROW constants should be defined globally or passed in
        const rawLength = Math.hypot(dx, dy);
        const scaledLength = rawLength / FORCE_ARROW.SCALE_RATIO;
        state.force.length = this.constrainValue(
            scaledLength,
            FORCE_ARROW.MIN_LENGTH,
            FORCE_ARROW.MAX_LENGTH
        );
    },

    updateNodePosition(pos) {
        // Ensure a valid node is selected (redundant check, but safe)
        if (state.selectedNode === null) return; // Removed the check for index 2 or 3

        const physicalPos = coordConversion.canvasToPhysical(pos, canvases.fem);
        const currentNode = state.nodes[state.selectedNode];

        // Apply constraints (e.g., keep within 0-1 physical space)
        currentNode.x = this.constrainValue(physicalPos.x, 0, 1);
        currentNode.y = this.constrainValue(physicalPos.y, 0, 1);
    },

    constrainValue(value, min, max) {
        return Math.max(min, Math.min(max, value));
    }
};

// Canvas management
const canvasManager = {
    calculateSize() {
        const windowHeight = window.innerHeight;
        if (windowHeight < 600) return 250;
        if (windowHeight < 800) return 300;
        //if (windowHeight < 1000) return 250;
        return 350;
    },

    resize(canvas) {
        if (!canvas) return;
        const size = this.calculateSize();
        canvas.width = size;
        canvas.height = size;
        canvas.style.width = `${size}px`;
        canvas.style.height = `${size}px`;
    },

    resizeAll() {
        this.resize(canvases.fem);
        this.resize(canvases.reference);
        this.resize(canvases.result);
        renderer.element();
        renderer.reference();
        const displacements = calculateDisplacements(state.nodes, state.force);
        renderer.results(displacements);
    },

    initialize() {
        // Get canvas elements
        canvases.fem = document.getElementById('femCanvas');
        canvases.reference = document.getElementById('refCanvas');
        canvases.result = document.getElementById('resultCanvas');

        if (!canvases.fem || !canvases.reference || !canvases.result) {
            console.error('Could not find all canvas elements');
            return;
        }

        // Get contexts
        canvases.contexts.fem = canvases.fem.getContext('2d');
        canvases.contexts.reference = canvases.reference.getContext('2d');
        canvases.contexts.result = canvases.result.getContext('2d');

        // Set initial sizes
        this.resizeAll(); // Call before adding listeners potentially

        // --- Add Event Listeners ---
        window.addEventListener('resize', () => this.resizeAll());

        // -- Mouse Listeners --
        canvases.fem.addEventListener('mousedown', (e) => handlers.FemMouseDown(e));
        canvases.fem.addEventListener('mousemove', (e) => handlers.FemMouseMove(e));
        canvases.fem.addEventListener('mouseup', (e) => handlers.FemMouseUp(e));
        canvases.fem.addEventListener('mouseleave', (e) => handlers.FemMouseUp(e)); // Treat leave as up

        canvases.reference.addEventListener('mousedown', (e) => handlers.RefMouseDown(e));
        canvases.reference.addEventListener('mousemove', (e) => handlers.RefMouseMove(e));
        canvases.reference.addEventListener('mouseup', (e) => handlers.RefMouseUp(e));
        canvases.reference.addEventListener('mouseleave', (e) => handlers.RefMouseUp(e)); // Treat leave as up

        document.addEventListener('mouseup', (e) => handlers.GlobalMouseUp(e)); // Global mouse up

        // -- Touch Listeners --
        // Use passive: false to allow preventDefault in move handlers
        canvases.fem.addEventListener('touchstart', (e) => handlers.FemMouseDown(e), { passive: true }); // Start doesn't need preventDefault immediately
        canvases.fem.addEventListener('touchmove', (e) => handlers.FemMouseMove(e), { passive: false });
        canvases.fem.addEventListener('touchend', (e) => handlers.FemMouseUp(e));
        canvases.fem.addEventListener('touchcancel', (e) => handlers.FemMouseUp(e)); // Treat cancel like end

        canvases.reference.addEventListener('touchstart', (e) => handlers.RefMouseDown(e), { passive: true }); // Start doesn't need preventDefault immediately
        canvases.reference.addEventListener('touchmove', (e) => handlers.RefMouseMove(e), { passive: false });
        canvases.reference.addEventListener('touchend', (e) => handlers.RefMouseUp(e));
        canvases.reference.addEventListener('touchcancel', (e) => handlers.RefMouseUp(e)); // Treat cancel like end

        document.addEventListener('touchend', (e) => handlers.GlobalMouseUp(e)); // Global touch end
        document.addEventListener('touchcancel', (e) => handlers.GlobalMouseUp(e)); // Global touch cancel


        // Initialize buttons and controls AFTER setting up canvas listeners
        controls.initialize();

        // Initial render
        renderer.element();
        renderer.reference();
        renderer.results(); // Initial calculation and rendering of the default view

        updateLessons(state.currentGaussPoint, state.nodes, state.force);
    }
};

window.addEventListener('load', () => {
    // Ensure initialization runs even if DOMContentLoaded fired early
    if (!canvases.fem) {
         console.log("Running initialization on window.load");
         canvasManager.initialize();
     }
});

document.addEventListener('DOMContentLoaded', () => {

    const toggleButton = document.getElementById('toggleControlsBtn');
    const controlsPanel = document.getElementById('controls'); // Needed for focus management potentially
    const body = document.body;

    // Check if elements exist
    if (!toggleButton || !controlsPanel) {
        console.warn("Toggle button or controls panel not found. Toggle functionality disabled.");
        return;
    }

    // Function to toggle the controls
    function toggleControls() {
        const isHidden = body.classList.contains('controls-hidden');

        body.classList.toggle('controls-hidden');

        // Update ARIA attribute and button label
        if (isHidden) {
            // Panel is now visible
            toggleButton.setAttribute('aria-expanded', 'true');
            toggleButton.setAttribute('aria-label', 'Hide Controls');
        } else {
            // Panel is now hidden
            toggleButton.setAttribute('aria-expanded', 'false');
            toggleButton.setAttribute('aria-label', 'Show Controls');
        }

        // Canvas Resizing
        setTimeout(() => {
             if (typeof canvasManager !== 'undefined' && canvasManager.resizeAll) {
                 canvasManager.resizeAll();
             }
        }, 350); // Match transition duration
    }

    // Set initial state based on default CSS (starts visible)
    toggleButton.setAttribute('aria-expanded', 'true');
    toggleButton.setAttribute('aria-label', 'Hide Controls');

    // Add click listener to the button
    toggleButton.addEventListener('click', toggleControls);

    // --- Initialize other things
    if (typeof canvasManager !== 'undefined' && canvasManager.initialize) {
       canvasManager.initialize();
    }
    // --- End Initialization ---

}); // End DOMContentLoaded