
// Enhancement Engine JavaScript Utilities

class EnhancementEngine {
    constructor() {
        this.initializeEventListeners();
    }
    
    initializeEventListeners() {
        // Form validation
        const forms = document.querySelectorAll('.needs-validation');
        forms.forEach(form => {
            form.addEventListener('submit', this.handleFormSubmit.bind(this));
        });
        
        // Gene example buttons
        document.addEventListener('click', (e) => {
            if (e.target.classList.contains('example-gene')) {
                this.handleGeneExampleClick(e);
            }
        });
    }
    
    handleFormSubmit(event) {
        const form = event.target;
        
        if (!form.checkValidity()) {
            event.preventDefault();
            event.stopPropagation();
        } else {
            this.showLoadingState(form);
        }
        
        form.classList.add('was-validated');
    }
    
    showLoadingState(form) {
        const submitBtn = form.querySelector('button[type="submit"]');
        if (submitBtn) {
            const originalText = submitBtn.innerHTML;
            submitBtn.innerHTML = '<span class="loading-spinner"></span> Analyzing...';
            submitBtn.disabled = true;
            
            // Store original text for potential restoration
            submitBtn.dataset.originalText = originalText;
        }
    }
    
    handleGeneExampleClick(event) {
        const button = event.target;
        const gene = button.dataset.gene;
        const variant = button.dataset.variant;
        
        if (gene) {
            const geneInput = document.getElementById('gene');
            const variantInput = document.getElementById('variant');
            
            if (geneInput) geneInput.value = gene;
            if (variantInput && variant) variantInput.value = variant;
            
            // Visual feedback
            this.highlightSelectedGene(button);
        }
    }
    
    highlightSelectedGene(selectedButton) {
        // Remove highlight from all gene buttons
        document.querySelectorAll('.example-gene').forEach(btn => {
            btn.classList.remove('selected');
        });
        
        // Highlight selected button
        selectedButton.classList.add('selected');
    }
    
    // Utility method for API calls
    async analyzeGene(geneData) {
        try {
            const response = await fetch('/api/analyze', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(geneData)
            });
            
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            
            const result = await response.json();
            return result;
        } catch (error) {
            console.error('Analysis failed:', error);
            throw error;
        }
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new EnhancementEngine();
});

// Export for use in other scripts
window.EnhancementEngine = EnhancementEngine;
