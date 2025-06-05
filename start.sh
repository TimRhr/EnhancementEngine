#!/bin/bash

# Enhancement Engine Docker Starter Script
# This script sets up and starts the Enhancement Engine with Docker

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ASCII Art Banner
echo -e "${BLUE}"
cat << "EOF"
╔═══════════════════════════════════════════════════════════════╗
║  🧬 Enhancement Engine - Docker Setup                         ║  
║                                                               ║
║  Comprehensive genetic enhancement simulation and analysis    ║
╚═══════════════════════════════════════════════════════════════╝
EOF
echo -e "${NC}"

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if Docker is installed
check_docker() {
    print_status "Checking Docker installation..."
    
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install Docker first:"
        echo "  - macOS: https://docs.docker.com/docker-for-mac/install/"
        echo "  - Windows: https://docs.docker.com/docker-for-windows/install/"
        echo "  - Linux: https://docs.docker.com/engine/install/"
        exit 1
    fi
    
    if ! command -v docker compose &> /dev/null; then
        print_error "Docker Compose is not installed. Please install Docker Compose:"
        echo "  - https://docs.docker.com/compose/install/"
        exit 1
    fi
    
    print_status "✅ Docker and Docker Compose are installed"
}

# Check if Docker daemon is running
check_docker_daemon() {
    print_status "Checking Docker daemon..."
    
    if ! docker info &> /dev/null; then
        print_error "Docker daemon is not running. Please start Docker first."
        exit 1
    fi
    
    print_status "✅ Docker daemon is running"
}

# Setup environment file
setup_env() {
    print_status "Setting up environment configuration..."
    
    if [ ! -f .env ]; then
        if [ -f .env.example ]; then
            cp .env.example .env
            print_status "✅ Created .env file from .env.example"
            print_warning "Please review and modify .env file if needed"
        else
            print_warning "No .env.example found, using default environment"
        fi
    else
        print_status "✅ .env file already exists"
    fi
}

# Build Docker images
build_images() {
    print_status "Building Docker images..."
    
    if docker compose build; then
        print_status "✅ Docker images built successfully"
    else
        print_error "Failed to build Docker images"
        exit 1
    fi
}

# Start services
start_services() {
    print_status "Starting Enhancement Engine services..."
    
    if docker compose up -d; then
        print_status "✅ Services started successfully"
    else
        print_error "Failed to start services"
        exit 1
    fi
}

# Wait for services to be ready
wait_for_services() {
    print_status "Waiting for services to be ready..."
    
    # Wait for the main app to be ready
    local max_attempts=30
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        if curl -f http://localhost:5000 &> /dev/null; then
            print_status "✅ Enhancement Engine is ready!"
            break
        fi
        
        if [ $attempt -eq $max_attempts ]; then
            print_error "Services did not start properly. Check logs with: docker compose logs"
            exit 1
        fi
        
        echo -n "."
        sleep 2
        ((attempt++))
    done
    echo ""
}

# Show service status
show_status() {
    print_status "Service Status:"
    docker compose ps
    echo ""
    
    print_status "Service URLs:"
    echo "  🌐 Web Interface: http://localhost:5000"
    echo "  🔧 Redis Cache:   http://localhost:6379 (if enabled)"
    echo ""
    
    print_status "Useful Commands:"
    echo "  📋 View logs:     docker compose logs -f"
    echo "  🛑 Stop services: docker compose down"
    echo "  🔄 Restart:       docker compose restart"
    echo "  🧹 Clean up:      docker compose down -v"
    echo ""
    
    print_status "Using Makefile:"
    echo "  📋 View logs:     make logs"
    echo "  🛑 Stop services: make down"
    echo "  🔄 Restart:       make restart"
    echo "  🧹 Clean up:      make clean"
    echo "  🆘 Help:          make help"
}

# Main execution
main() {
    echo -e "${GREEN}Starting Enhancement Engine setup...${NC}\n"
    
    # Perform all checks and setup
    check_docker
    check_docker_daemon
    setup_env
    build_images
    start_services
    wait_for_services
    
    # Show final status
    echo -e "\n${GREEN}🎉 Enhancement Engine is now running!${NC}\n"
    show_status
    
    echo -e "${BLUE}Happy gene editing! 🧬${NC}"
}

# Handle script arguments
case "${1:-}" in
    "build")
        check_docker
        check_docker_daemon
        build_images
        ;;
    "start")
        check_docker
        check_docker_daemon
        start_services
        wait_for_services
        show_status
        ;;
    "stop")
        print_status "Stopping Enhancement Engine..."
        docker compose down
        print_status "✅ Services stopped"
        ;;
    "status")
        show_status
        ;;
    "clean")
        print_status "Cleaning up Enhancement Engine..."
        docker compose down -v --remove-orphans
        print_status "✅ Cleanup complete"
        ;;
    "help"|"-h"|"--help")
        echo "Enhancement Engine Docker Starter"
        echo ""
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  (no args)  - Full setup and start"
        echo "  build      - Only build images"
        echo "  start      - Only start services"
        echo "  stop       - Stop all services"
        echo "  status     - Show service status"
        echo "  clean      - Stop and remove all data"
        echo "  help       - Show this help"
        ;;
    *)
        main
        ;;
esac