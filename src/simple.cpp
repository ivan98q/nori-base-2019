#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
protected:
  Point3f m_position;
  Color3f m_energy;

public:
  SimpleIntegrator(const PropertyList &props) {
    m_position = props.getPoint("position");
    m_energy = props.getColor("energy");
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
    Intersection its;
    // If it doesn't hit the scene don't render it 
    if (!scene->rayIntersect(ray,its))
      return Color3f(0.0f);

    Vector3f direction_x_p = m_position - its.p;
    Ray3f visibleRay(its.p,direction_x_p);
    if(scene->rayIntersect(visibleRay)){
      return Color3f(0.0f);
    }
    Color3f color = Color3f(1.0f);

    float theta = its.shFrame.toLocal(direction_x_p).normalized()[2]; // This is the theta that's the normal between x and p. We only need z.
    // m_PI comes from math but it's already defined somewhere!
    color = (m_energy / (4.0f * (M_PI*M_PI))) *(std::max(0.0f,theta)/direction_x_p.dot(direction_x_p)) * color; // This is the equation
    return color;
  }

  std::string toString() const {
    return "SimpleIntegrator[]";
  }
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END
