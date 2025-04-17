export class MatrixFormatter {
    constructor(options = {}) {
        this.defaultOptions = {
            precision: 4,
            minWidth: 8,
            useSymbolicMinus: true,
            className: 'matrix'
        };
        this.options = { ...this.defaultOptions, ...options };
    }

    formatNumber(num, { isStiffness = false, precision } = {}) {
        if (typeof num !== 'number') return num;

        if (isStiffness) {
            const intValue = Math.round(num);
            const numStr = this.options.useSymbolicMinus ?
                intValue.toString().replace('-', '−') :
                intValue.toString();
            return numStr.padStart(this.options.minWidth + (intValue < 0 ? 0 : 1), ' ');
        }

        const usePrecision = precision !== undefined ? precision : this.options.precision;
        const formatted = num.toFixed(usePrecision);
        return this.options.useSymbolicMinus ?
            formatted.replace(/^-/, '−') :
            formatted;
    }

    createMatrix(data, config = {}) {
        const {
            colHeaders = [],
            rowHeaders = [],
            title = '',
            isStiffness = false,
            customClasses = [],
            precision
        } = config;

        const classes = [this.options.className, ...customClasses].join(' ');
        let html = '<div class="matrix-container">';

        if (title) {
            html += `<h5>${title}</h5>`;
        }

        html += `<div class="${classes}">`;

        // Add column headers if present
        if (colHeaders.length > 0) {
            html += '<div class="matrix-row">';
            // Add empty cell for row header column if row headers exist
            if (rowHeaders.length > 0) {
                html += '<div class="matrix-cell matrix-header"></div>';
            }
            colHeaders.forEach(header => {
                const formattedHeader = isStiffness ?
                    header.padStart(this.options.minWidth, ' ') :
                    header;
                html += `<div class="matrix-cell matrix-header">${formattedHeader}</div>`;
            });
            html += '</div>';
        }

        // Process matrix data
        const matrix = Array.isArray(data[0]) ? data : [data];
        matrix.forEach((row, i) => {
            html += '<div class="matrix-row">';

            // Add row header if present
            if (rowHeaders.length > i) {
                const formattedHeader = isStiffness ?
                    rowHeaders[i].padStart(this.options.minWidth, ' ') :
                    rowHeaders[i];
                html += `<div class="matrix-cell matrix-header">${formattedHeader}</div>`;
            }

            // Add data cells
            row.forEach(cell => {
                html += `<div class="matrix-cell">${this.formatNumber(cell, {
                    isStiffness: isStiffness,
                    precision: precision
                })}</div>`;
            });

            html += '</div>';
        });

        html += '</div></div>';
        return html;
    }

    formatGaussPoints(points, weights, precision) {
        const data = points.map((point, i) => [...point, weights[i]]);
        return this.createMatrix(data, {
            colHeaders: ['ξ', 'η', 'Weight'],
            rowHeaders: points.map((_, i) => `Point ${i + 1}`),
            title: 'Gauss Points and Weights',
            precision
        });
    }
}

// Predefined matrix configurations
export const matrixConfigs = {
    'N': {
        colHeaders: ['N₁', 'N₂', 'N₃', 'N₄'],
        rowHeaders: ['Shape Functions']
    },
    'dN/dξ': {
        colHeaders: ['∂N/∂ξ', '∂N/∂η'],
        rowHeaders: ['N₁', 'N₂', 'N₃', 'N₄']
    },
    'J': {
        colHeaders: ['∂x', '∂y'],
        rowHeaders: ['/∂ξ', '/∂η']
    },
     'J⁻¹': {
        colHeaders: ['', ''],
        rowHeaders: ['', '']
    },
    'dN/dx': {
        colHeaders: ['∂N/∂x', '∂N/∂y'],
        rowHeaders: ['N₁', 'N₂', 'N₃', 'N₄']
    },
    'B': {
        colHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
        rowHeaders: ['εxx', 'εyy', 'εxy']
    },
    'D': {
        colHeaders: ['σxx', 'σyy', 'σxy'],
        rowHeaders: ['εxx', 'εyy', 'εxy']
    },
    'K': {
        colHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
        rowHeaders: ['u₁', 'v₁', 'u₂', 'v₂', 'u₃', 'v₃', 'u₄', 'v₄'],
        isStiffness: true
    }
};



// Create default formatter instance
export const formatter = new MatrixFormatter();