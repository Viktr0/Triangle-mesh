// -*- mode: c++ -*-
#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

using qglviewer::Vec;

class MyViewer : public QGLViewer {
  Q_OBJECT

public:
  explicit MyViewer(QWidget *parent);
  virtual ~MyViewer();

  inline double getCutoffRatio() const;
  inline double X_getMedian() const;
  inline void setCutoffRatio(double ratio);
  inline double getMeanMin() const;
  inline void setMeanMin(double min);
  inline double getMeanMax() const;
  inline void setMeanMax(double max);
  inline const double *getSlicingDir() const;
  inline void setSlicingDir(double x, double y, double z);
  inline double getSlicingScaling() const;
  inline void setSlicingScaling(double scaling);
  bool openMesh(const std::string &filename, bool update_view = true);
  bool openBezier(const std::string &filename, bool update_view = true);
  bool saveBezier(const std::string &filename);
  bool openCoons(const std::string &filename, bool update_view = true);
  bool saveCoons(const std::string &filename);

signals:
  void startComputation(QString message);
  void midComputation(int percent);
  void endComputation();
  void X_showMedian();

protected:
  virtual void init() override;
  virtual void draw() override;
  virtual void drawWithNames() override;
  virtual void postSelection(const QPoint &p) override;
  virtual void keyPressEvent(QKeyEvent *e) override;
  virtual void mouseMoveEvent(QMouseEvent *e) override;
  virtual QString helpString() const override;

private:
  struct MyTraits : public OpenMesh::DefaultTraits {
    using Point  = OpenMesh::Vec3d; // the default would be Vec3f
    using Normal = OpenMesh::Vec3d;
    VertexTraits {
		double mean;              // approximated mean curvature
		double u;
		double v;
    };
  };
  using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
  using Vector = OpenMesh::VectorT<double,3>;

  // Mesh
  void updateMesh(bool update_mean_range = true);
  void updateVertexNormals();
#ifdef USE_JET_FITTING
  void updateWithJetFit(size_t neighbors);
#endif
  void localSystem(const Vector &normal, Vector &u, Vector &v);
  double voronoiWeight(MyMesh::HalfedgeHandle in_he);
  void updateMeanMinMax();
  void updateMeanCurvature();
  void updateContinuousMeanCurvature();

  // Bezier
  static void bernsteinAll(size_t n, double u, std::vector<double> &coeff);
  void generateMesh(size_t resolution);

  // Coons
  size_t findSpanA(double u, int edge_idx) const;
  void basisFunctionsA(size_t i, double u, std::vector<double> &coeff, int edge_idx, int p) const;
  Vec evaluateA(double u, int edge_idx, int p, std::vector<Vec> cps) const;
  void derivativeControlPoints(size_t d, size_t r1, size_t r2, std::vector<std::vector<Vec>> &dcp, std::vector<Vec> cp, int edge_idx) const;
  void generateMeshCoons(size_t resolution);
  std::vector<Vec> getControlPoints(int edge_idx) const;


  // BSpline
  struct BSplineCurve
  {
	  size_t p = 3;						// degree
	  size_t n;						// n + 1 = cp.size()
	  std::vector<double> knots;      // first and last p+1 values are the same ("clamped")
	  std::vector<Vec> cp;			// knots.size() = cp.size() + p + 1

	  void updateCP() const;
	  size_t findSpan(double u) const;
	  void basisFunctions(size_t i, double u, std::vector<double> &coeff) const;
	  Vec evaluate(double u) const;
	  void basisFunctionDerivatives(size_t i, double u, size_t d, std::vector<std::vector<double>> &der) const;
	  Vec derivatives(double u, size_t d, std::vector<Vec> &der) const;
	  void derivativeControlPoints(size_t d, size_t r1, size_t r2, std::vector<std::vector<Vec>> &dcp) const;
	  void basisFunctionsAll(size_t i, double u, std::vector<std::vector<double>> &coeff) const;
	  Vec derivativesByControlPoints(double u, size_t d, std::vector<Vec> &der) const;
  };

  // Visualization
  void setupCamera();
  Vec meanMapColor(double d) const;
  void drawControlNet() const;
  void drawControlNetCoons() const;
  void drawAxes() const;
  void drawAxesWithNames() const;
  static Vec intersectLines(const Vec &ap, const Vec &ad, const Vec &bp, const Vec &bd);
  
  // Task
  double X_getTriangleArea(OpenMesh::FaceHandle face) const;
  std::vector<double> X_triangleAreas() const;
  double X_calcMedian(std::vector<double> trianglesArray) const;
  double X_calcFirstQuartile(std::vector<double> trianglesArray) const;
  double X_calcLastQuartile(std::vector<double> trianglesArray) const;
  std::vector<GLdouble> X_triangleColor(OpenMesh::FaceHandle face, double firstQuartile, double lastQuartile) const;

  // Other
  void fairMesh();

  //////////////////////
  // Member variables //
  //////////////////////

  enum class ModelType { NONE, MESH, BEZIER_SURFACE, COONS_SURFACE } model_type;

  // Mesh
  MyMesh mesh;

  // Bezier
  size_t degree[2];
  std::vector<Vec> control_points;

  // Coons
  BSplineCurve bs[4];				// oldalak
  size_t edge_info[5];			// oldalak kezdo indexe es pontok szama
  size_t p;						// fokszam
  size_t segments[4];				// n + 1 = control_p.size()
  std::vector<double> knots[4];	// elso es utolso p+1 ertek megegyezik
  std::vector<Vec> testBSpoint;
  //std::vector<Vec> control_p;	// knots.size() = control_p.size() + p + 1

  // Visualization
  double mean_min, mean_max, cutoff_ratio;
  double X_firstQuartile, X_lastQuartile;
  bool show_control_points, show_solid, show_wireframe;
  enum class Visualization { PLAIN, MEAN, SLICING, ISOPHOTES, X_TRIANGLES, CONTINUOUS_MEAN } visualization;
  GLuint isophote_texture, environment_texture, current_isophote_texture, slicing_texture;
  Vector slicing_dir;
  double slicing_scaling;
  int selected_vertex;
  struct ModificationAxes {
    bool shown;
    float size;
    int selected_axis;
    Vec position, grabbed_pos, original_pos;
  } axes;
  std::string last_filename;
};

#include "MyViewer.hpp"
