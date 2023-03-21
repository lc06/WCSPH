#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <Eigen/Dense>
#include <EigenExp/CXX11/Tensor>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <learnopengl//shader.h>
#include <config.h>
#include <particle.h>
#include <iostream>
#include <algorithm>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void initParticles();
void storeToVertices();
std::vector<int> findNeighbors(Particle* p);
void clearForce();
void step();
void updateDt();

glm::ivec2 resolution = glm::ivec2(800, 800);

std::vector<Particle> boundaryParticles(boundaryParticleNum, Particle());
std::vector<Particle> fluidParticles(fluidParticleNum, Particle());
std::vector<Particle> obstacleParticles(obstacleParticleNum, Particle());

std::vector<float>vertices (4 * (boundaryParticleNum + fluidParticleNum + obstacleParticleNum), 0.0);

int grid_dim = int(12.7 / grid_size) + 1;
std::vector<std::vector<int>> grids (grid_dim * grid_dim, std::vector<int> {});
int offset[][2] = {{-1, -1}, {0, -1}, {1, -1},
                   {-1, 0}, {0, 0}, {1, 0},
                   {-1, 1}, {0, 1}, {1, 1}};

bool isFirst = true;
int frame = 0;

int main() {

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(resolution[0], resolution[1], "SPH", nullptr, nullptr);
    if(window == nullptr) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    if(!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    Shader myShader((path + "resources\\shaders\\sph.vert").c_str(), (path + "resources\\shaders\\sph.frag").c_str());

    initParticles();
    storeToVertices();

    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STREAM_DRAW);

//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
//    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*) nullptr);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2);

    myShader.use();
    myShader.setVec2("u_resolution", resolution);

    glm::mat4 model(1.0);
    model = glm::translate(model, glm::vec3(-0.89, -0.89, 0.0));
    model = glm::scale(model, glm::vec3(1. / 7., 1. / 7., 1.));
    myShader.setMat4("model", model);

    while(!glfwWindowShouldClose(window)) {

        processInput(window);

        glClear(GL_COLOR_BUFFER_BIT);
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

        step();
        updateDt();
        storeToVertices();

        std::cout << frame << std::endl;
        frame++;

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STREAM_DRAW);

        myShader.use();
        float timeValue = glfwGetTime();
        myShader.setFloat("u_time", timeValue);

        glPointSize(2);
        glBindVertexArray(VAO);
        glDrawArrays(GL_POINTS, 0, boundaryParticleNum + fluidParticleNum + obstacleParticleNum);
//        glDrawElements(GL_POINTS, 6, GL_UNSIGNED_INT, nullptr);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);

    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
    resolution[0] = width;
    resolution[1] = height;
}

void processInput(GLFWwindow* window) {
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void initParticles() {

    for (int i = 0; i < boundaryParticleNum; i+=4) {
        boundaryParticles[i].set_p_x(0.1 + i / 4 * padding);
        boundaryParticles[i].set_p_y(0.1);
        boundaryParticles[i].set_type(0);

        boundaryParticles[i + 1].set_p_x(12.7);
        boundaryParticles[i + 1].set_p_y(0.1 + i / 4 * padding);
        boundaryParticles[i + 1].set_type(0);

        boundaryParticles[i + 2].set_p_x(12.7 - i / 4 * padding);
        boundaryParticles[i + 2].set_p_y(12.7);
        boundaryParticles[i + 2].set_type(0);

        boundaryParticles[i + 3].set_p_x(0.1);
        boundaryParticles[i + 3].set_p_y(12.7 - i / 4 * padding);
        boundaryParticles[i + 3].set_type(0);
    }

    // 最上层实际有31个粒子
    int w = 60, cnt = 0, y=0;
    for(int j = 0; j < fluidParticleNum; j++) {
        fluidParticles[j].set_p_x(cnt % w * padding + x_bias);
        fluidParticles[j].set_p_y(y * padding + y_bias);
        fluidParticles[j].set_type(1);
        fluidParticles[j].set_id(j);

        cnt ++;
        if(cnt % w == 0)
            y++;
        if(cnt % (2 * w) == 0) {
            w -= 1;
            cnt = 0;
        }
    }

    for(int k=0; k < obstacleParticleNum; k++) {
        double degree = 3 * k * pi / 180.0;
        double p_x = cos(degree) + 10.0;
        double p_y = sin(degree) + 0.1;
        obstacleParticles[k].set_p_x(p_x);
        obstacleParticles[k].set_p_y(p_y);
        obstacleParticles[k].set_type(2);
    }
}

void storeToVertices() {
    int i =0, j=0, k=0;
    if(isFirst)
        for(; i< boundaryParticleNum; i++) {
            vertices[i*4] = boundaryParticles[i].get_p_x();
            vertices[i*4+1] = boundaryParticles[i].get_p_y();
            vertices[i*4+2] = boundaryParticles[i].get_v();
            vertices[i*4+3] = boundaryParticles[i].get_type();
        }
    else
        i = boundaryParticleNum;

    for(; j< fluidParticleNum; j++) {
        vertices[(i+j)*4] = fluidParticles[j].get_p_x();
        vertices[(i+j)*4+1] = fluidParticles[j].get_p_y();
        vertices[(i+j)*4+2] = fluidParticles[j].get_v();
        vertices[(i+j)*4+3] = fluidParticles[j].get_type();
    }

    if(isFirst)
        for(; k< obstacleParticleNum; k++) {
            vertices[(i+j+k)*4] = obstacleParticles[k].get_p_x();
            vertices[(i+j+k)*4+1] = obstacleParticles[k].get_p_y();
            vertices[(i+j+k)*4+2] = obstacleParticles[k].get_v();
            vertices[(i+j+k)*4+3] = obstacleParticles[k].get_type();
        }
    else isFirst = false;
}

void calcGrids() {

//    grids.clear();

    for(auto& subVec: grids) subVec.clear();

    int i =0, j=0, k=0;
    for (; i < boundaryParticleNum; i++) {
        int grid_x = int(boundaryParticles[i].get_p_x() / grid_size);
        int grid_y = int(boundaryParticles[i].get_p_y() / grid_size);

        grids[grid_y * grid_dim + grid_x].push_back(i);
    }

    for (; j < fluidParticleNum; j++) {
        int grid_x = int(fluidParticles[j].get_p_x() / grid_size);
        int grid_y = int(fluidParticles[j].get_p_y() / grid_size);

        grids[grid_y * grid_dim + grid_x].push_back(i + j);
    }

    for (; k < obstacleParticleNum; k++) {
        int grid_x = int(obstacleParticles[k].get_p_x() / grid_size);
        int grid_y = int(obstacleParticles[k].get_p_y() / grid_size);

        grids[grid_y * grid_dim + grid_x].push_back(i + j + k);
    }
}

std::vector<int> findNeighbors(Particle* p) {
    int grid_x = int(p->get_p_x() / grid_size);
    int grid_y = int(p->get_p_y() / grid_size);

    std::vector<int> res {};
    for(auto& dire: offset) {
        int offset_grid_x = grid_x + dire[0];
        int offset_grid_y = grid_y + dire[1];

        if(offset_grid_x < 0 || offset_grid_x >= grid_dim || offset_grid_y < 0 || offset_grid_y >= grid_dim) continue;

        auto grid_dire = grids[offset_grid_y * grid_dim + offset_grid_x];
        res.insert(res.end(), grid_dire.begin(), grid_dire.end());
    }
    return res;
}

void clearForce() {
    for (auto &p : fluidParticles) {
        p.clear_a_x();
        p.clear_a_y();
    }
}

void step() {

    calcGrids();
    clearForce();

    for (auto &p : fluidParticles) {
        std::vector<int> neighborsIndex = findNeighbors(&p);
        std::vector<Particle> neighbors;

        for (auto &index : neighborsIndex) {
            if(index < boundaryParticleNum) neighbors.push_back(boundaryParticles[index]);
            else if(index < boundaryParticleNum + fluidParticleNum)
              neighbors.push_back(fluidParticles[index - boundaryParticleNum]);
            else neighbors.push_back(obstacleParticles[index - boundaryParticleNum - fluidParticleNum]);
        }

        for(int j=0; j < neighbors.size(); j++) {
            if(neighbors[j].get_id() == p.get_id())
              continue;
            if(neighbors[j].get_type() == 1)
                p.update_rho(&neighbors[j]);
        }

        p.update_pressure();
        p.calc_gravity_force_acce();

        for(int j=0; j < neighbors.size(); j++) {

            if(neighbors[j].get_id() == p.get_id()) continue;
            if(neighbors[j].get_type() == 1) {
                p.calc_pressure_force_acce(&neighbors[j]);
                p.calc_viscosity_force_acce(&neighbors[j]);
            } else {
                p.calc_repulsive_force_acce(&neighbors[j]);
            }
        }

//        p.calc_acceleration();
        p.calc_velocity();
        p.calc_position();

        for(int j=0; j < neighbors.size(); j++) {
            if(neighbors[j].get_id() == p.get_id()) continue;
            if(neighbors[j].get_type() == 1)
                p.calc_xsph_correction(&neighbors[j]);
        }
    }
}

void updateDt() {
    double max_vel = 0.0, max_acce = 0.0, max_rho = 0.0;

    for (auto &p : fluidParticles) {
        if(p.get_v() > max_vel) max_vel = p.get_v();
        if(p.get_a() > max_acce) max_acce = p.get_a();
        if(p.get_rho() > max_rho) max_rho = p.get_rho();
    }

    std::cout << "dt: " << dt << ", max_vel: " << max_vel << std::endl;

    double dtCfl = CFL_V * h / max_vel;
    double dtF = CFL_A * sqrt(h / max_acce);
    double dtA = 0.20 * h / (c_0 * sqrt(pow((max_rho / rho_0), gamma_)));

    dt = dtCfl < dtF ? (dtCfl < dtA ? dtCfl : dtA) : (dtF < dtA ? dtF : dtA);
}