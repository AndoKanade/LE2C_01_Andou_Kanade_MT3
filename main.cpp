#define NOMINMAX
#include <Novice.h>
#include <windows.h>
#define _USE_MATH_DEFINES
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <imgui.h>
#include <math.h>

const char kWindowTitle[] = "LE2C_01_アンドウ_カナデ_MT3";

#pragma region 構造体
typedef struct Matrix4x4 {

  float m[4][4];
} Matrix4x4;

typedef struct Vector3 {
  float x;
  float y;
  float z;

  // 内積
  float Dot(const Vector3 &other) const {
    return x * other.x + y * other.y + z * other.z;
  }

  Vector3 operator-(const Vector3 &rhs) const {
    return {x - rhs.x, y - rhs.y, z - rhs.z};
  }

  Vector3 operator+(const Vector3 &rhs) const {
    return {x + rhs.x, y + rhs.y, z + rhs.z};
  }

  // スカラー倍
  Vector3 operator*(float scalar) const {
    return {x * scalar, y * scalar, z * scalar};
  }

  // ベクトルの長さの2乗（省略できる）
  float LengthSquared() const { return x * x + y * y + z * z; }

} Vector3;

typedef struct Sphere {
  Vector3 center;
  float radius;
  unsigned int color; // 色を追加
} Sphere;

struct Line {
  Vector3 origin;
  Vector3 diff;
};

struct Ray {
  Vector3 origin;
  Vector3 diff;
};

struct Segment {
  Vector3 origin;
  Vector3 diff;
};

struct Plane {
  Vector3 normal;
  float distance;
};

struct Camera {
  Vector3 position;
  float yaw;   // 左右
  float pitch; // 上下
};

#pragma endregion

#pragma region 関数

float Dot(const Vector3 &a, const Vector3 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 Add(const Vector3 &v1, const Vector3 &v2) {
  Vector3 result;

  result.x = v1.x + v2.x;
  result.y = v1.y + v2.y;
  result.z = v1.z + v2.z;

  return result;
}

Vector3 Subtract(const Vector3 &v1, const Vector3 &v2) {
  Vector3 result;

  result.x = v1.x - v2.x;
  result.y = v1.y - v2.y;
  result.z = v1.z - v2.z;

  return result;
}

Vector3 Multiply(float scalar, const Vector3 &v) {
  Vector3 result;

  result.x = scalar * v.x;
  result.y = scalar * v.y;
  result.z = scalar * v.z;

  return result;
}

Vector3 Normalize(const Vector3 &v) {
  float length = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  if (length == 0.0f)
    return {0, 0, 0};
  return {v.x / length, v.y / length, v.z / length};
}

Matrix4x4 Add(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;
  result.m[0][0] = m1.m[0][0] + m2.m[0][0];
  result.m[0][1] = m1.m[0][1] + m2.m[0][1];
  result.m[0][2] = m1.m[0][2] + m2.m[0][2];
  result.m[0][3] = m1.m[0][3] + m2.m[0][3];
  result.m[1][0] = m1.m[1][0] + m2.m[1][0];
  result.m[1][1] = m1.m[1][1] + m2.m[1][1];
  result.m[1][2] = m1.m[1][2] + m2.m[1][2];
  result.m[1][3] = m1.m[1][3] + m2.m[1][3];
  result.m[2][0] = m1.m[2][0] + m2.m[2][0];
  result.m[2][1] = m1.m[2][1] + m2.m[2][1];
  result.m[2][2] = m1.m[2][2] + m2.m[2][2];
  result.m[2][3] = m1.m[2][3] + m2.m[2][3];
  result.m[3][0] = m1.m[3][0] + m2.m[3][0];
  result.m[3][1] = m1.m[3][1] + m2.m[3][1];
  result.m[3][2] = m1.m[3][2] + m2.m[3][2];
  result.m[3][3] = m1.m[3][3] + m2.m[3][3];
  return result;
}

Matrix4x4 Subtract(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;
  result.m[0][0] = m1.m[0][0] - m2.m[0][0];
  result.m[0][1] = m1.m[0][1] - m2.m[0][1];
  result.m[0][2] = m1.m[0][2] - m2.m[0][2];
  result.m[0][3] = m1.m[0][3] - m2.m[0][3];
  result.m[1][0] = m1.m[1][0] - m2.m[1][0];
  result.m[1][1] = m1.m[1][1] - m2.m[1][1];
  result.m[1][2] = m1.m[1][2] - m2.m[1][2];
  result.m[1][3] = m1.m[1][3] - m2.m[1][3];
  result.m[2][0] = m1.m[2][0] - m2.m[2][0];
  result.m[2][1] = m1.m[2][1] - m2.m[2][1];
  result.m[2][2] = m1.m[2][2] - m2.m[2][2];
  result.m[2][3] = m1.m[2][3] - m2.m[2][3];
  result.m[3][0] = m1.m[3][0] - m2.m[3][0];
  result.m[3][1] = m1.m[3][1] - m2.m[3][1];
  result.m[3][2] = m1.m[3][2] - m2.m[3][2];
  result.m[3][3] = m1.m[3][3] - m2.m[3][3];
  return result;
}

Matrix4x4 Multiply(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result{};
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 4; ++col) {
      result.m[row][col] = 0.0f;
      for (int k = 0; k < 4; ++k) {
        result.m[row][col] += m1.m[row][k] * m2.m[k][col];
      }
    }
  }
  return result;
}

Matrix4x4 MakeTranslateMatrix(const Vector3 &translate) {
  Matrix4x4 matrix = {}; // すべて0で初期化
  // 単位行列の形に設定
  matrix.m[0][0] = 1.0f;
  matrix.m[1][1] = 1.0f;
  matrix.m[2][2] = 1.0f;
  matrix.m[3][3] = 1.0f;
  // 平行移動成分を設定
  matrix.m[3][0] = translate.x;
  matrix.m[3][1] = translate.y;
  matrix.m[3][2] = translate.z;
  return matrix;
}

Matrix4x4 MakeScaleMatrix(const Vector3 &scale) {
  Matrix4x4 matrix = {}; // すべて0で初期化
  // スケール行列の設定
  matrix.m[0][0] = scale.x;
  matrix.m[1][1] = scale.y;
  matrix.m[2][2] = scale.z;
  matrix.m[3][3] = 1.0f;
  return matrix;
}

Matrix4x4 MakeRotateXMatrix(float radian) {
  Matrix4x4 result{};

  result.m[0][0] = 1;
  result.m[3][3] = 1;

  // X軸回転に必要な部分だけ上書き
  result.m[1][1] = std::cos(radian);
  result.m[1][2] = std::sin(radian);
  result.m[2][1] = -std::sin(radian);
  result.m[2][2] = std::cos(radian);

  return result;
}

Matrix4x4 MakeRotateYMatrix(float radian) {
  Matrix4x4 result{};

  result.m[1][1] = 1.0f;
  result.m[3][3] = 1.0f;

  result.m[0][0] = std::cos(radian);
  result.m[0][2] = -std::sin(radian);
  result.m[2][0] = std::sin(radian);
  result.m[2][2] = std::cos(radian);

  return result;
}

Matrix4x4 MakeRotateZMatrix(float radian) {
  Matrix4x4 result{};

  result.m[2][2] = 1;
  result.m[3][3] = 1;

  result.m[0][0] = std::cos(radian);
  result.m[0][1] = std::sin(radian);
  result.m[1][0] = -std::sin(radian);
  result.m[1][1] = std::cos(radian);

  return result;
}

Matrix4x4 MakeAffineMatrix(const Vector3 &scale, const Vector3 &rotate,
                           const Vector3 &translate) {

  Matrix4x4 scaleMatrix = MakeScaleMatrix(scale);
  Matrix4x4 rotateX = MakeRotateXMatrix(rotate.x);
  Matrix4x4 rotateY = MakeRotateYMatrix(rotate.y);
  Matrix4x4 rotateZ = MakeRotateZMatrix(rotate.z);

  // 回転順: Z → X → Y →（スケーリング）→ 平行移動
  Matrix4x4 rotateMatrix = Multiply(Multiply(rotateX, rotateY), rotateZ);

  Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);

  Matrix4x4 affineMatrix =
      Multiply(Multiply(scaleMatrix, rotateMatrix), translateMatrix);

  return affineMatrix;
}

Vector3 Transform(const Vector3 &vector, const Matrix4x4 &matrix) {
  Vector3 result;
  result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] +
             vector.z * matrix.m[2][0] + matrix.m[3][0];
  result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] +
             vector.z * matrix.m[2][1] + matrix.m[3][1];
  result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] +
             vector.z * matrix.m[2][2] + matrix.m[3][2];
  float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] +
            vector.z * matrix.m[2][3] + matrix.m[3][3];
  assert(w != 0.0f);
  result.x /= w;
  result.y /= w;
  result.z /= w;

  return result;
}

Matrix4x4 Inverse(const Matrix4x4 &m) {
  Matrix4x4 result;

  // 行列の行列式を計算
  float det =
      m.m[0][0] *
          (m.m[1][1] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
           m.m[1][2] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) +
           m.m[1][3] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1])) -
      m.m[0][1] *
          (m.m[1][0] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
           m.m[1][2] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
           m.m[1][3] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0])) +
      m.m[0][2] *
          (m.m[1][0] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) -
           m.m[1][1] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
           m.m[1][3] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) -
      m.m[0][3] * (m.m[1][0] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1]) -
                   m.m[1][1] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0]) +
                   m.m[1][2] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0]));

  if (det == 0) {
    // 行列式がゼロの場合、逆行列は存在しません
    return result; // 逆行列は存在しないのでゼロ行列を返す
  }

  // 行列式の逆数を計算
  float invDet = 1.0f / det;

  // 各要素を余因子行列から計算
  result.m[0][0] =
      (m.m[1][1] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
       m.m[1][2] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) +
       m.m[1][3] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1])) *
      invDet;
  result.m[0][1] =
      (-m.m[0][1] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) +
       m.m[0][2] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) -
       m.m[0][3] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1])) *
      invDet;
  result.m[0][2] =
      (m.m[0][1] * (m.m[1][2] * m.m[3][3] - m.m[1][3] * m.m[3][2]) -
       m.m[0][2] * (m.m[1][1] * m.m[3][3] - m.m[1][3] * m.m[3][1]) +
       m.m[0][3] * (m.m[1][1] * m.m[3][2] - m.m[1][2] * m.m[3][1])) *
      invDet;
  result.m[0][3] =
      (-m.m[0][1] * (m.m[1][2] * m.m[2][3] - m.m[1][3] * m.m[2][2]) +
       m.m[0][2] * (m.m[1][1] * m.m[2][3] - m.m[1][3] * m.m[2][1]) -
       m.m[0][3] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1])) *
      invDet;

  result.m[1][0] =
      (-m.m[1][0] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) +
       m.m[1][2] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) -
       m.m[1][3] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0])) *
      invDet;
  result.m[1][1] =
      (m.m[0][0] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
       m.m[0][2] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
       m.m[0][3] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0])) *
      invDet;
  result.m[1][2] =
      -(m.m[0][0] * (m.m[1][2] * m.m[3][3] - m.m[1][3] * m.m[3][2]) -
        m.m[0][2] * (m.m[1][0] * m.m[3][3] - m.m[1][3] * m.m[3][0]) +
        m.m[0][3] * (m.m[1][0] * m.m[3][2] - m.m[1][2] * m.m[3][0])) *
      invDet;

  result.m[1][3] =
      (m.m[0][0] * (m.m[1][2] * m.m[2][3] - m.m[1][3] * m.m[2][2]) -
       m.m[0][2] * (m.m[1][0] * m.m[2][3] - m.m[1][3] * m.m[2][0]) +
       m.m[0][3] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0])) *
      invDet;

  result.m[2][0] =
      (m.m[1][0] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) -
       m.m[1][1] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
       m.m[1][3] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[2][1] =
      (-m.m[0][0] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) +
       m.m[0][1] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) -
       m.m[0][3] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[2][2] =
      (m.m[0][0] * (m.m[1][1] * m.m[3][3] - m.m[1][3] * m.m[3][1]) -
       m.m[0][1] * (m.m[1][0] * m.m[3][3] - m.m[1][3] * m.m[3][0]) +
       m.m[0][3] * (m.m[1][0] * m.m[3][1] - m.m[1][1] * m.m[3][0])) *
      invDet;

  result.m[2][3] =
      (-m.m[0][0] * (m.m[1][1] * m.m[2][3] - m.m[1][3] * m.m[2][1]) +
       m.m[0][1] * (m.m[1][0] * m.m[2][3] - m.m[1][3] * m.m[2][0]) -
       m.m[0][3] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0])) *
      invDet;

  result.m[3][0] =
      (-m.m[1][0] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1]) +
       m.m[1][1] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0]) -
       m.m[1][2] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[3][1] =
      (m.m[0][0] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1]) -
       m.m[0][1] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0]) +
       m.m[0][2] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[3][2] =
      (-m.m[0][0] * (m.m[1][1] * m.m[3][2] - m.m[1][2] * m.m[3][1]) +
       m.m[0][1] * (m.m[1][0] * m.m[3][2] - m.m[1][2] * m.m[3][0]) -
       m.m[0][2] * (m.m[1][0] * m.m[3][1] - m.m[1][1] * m.m[3][0])) *
      invDet;

  result.m[3][3] =
      (m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
       m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
       m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0])) *
      invDet;

  return result;
}

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio,
                                   float nearClip, float farClip) {

  float f = 1.0f / tanf(fovY * 0.5f);
  float range = farClip / (farClip - nearClip);

  Matrix4x4 result = {};

  result.m[0][0] = f / aspectRatio;
  result.m[1][1] = f;
  result.m[2][2] = range;
  result.m[2][3] = 1.0f;
  result.m[3][2] = -range * nearClip;

  return result;
}

Matrix4x4 MakeOrthographicMatrix(float left, float top, float right,
                                 float bottom, float nearClip, float farClip) {
  Matrix4x4 result = {};

  result.m[0][0] = 2.0f / (right - left);
  result.m[1][1] = 2.0f / (top - bottom);
  result.m[2][2] = 1.0f / (farClip - nearClip);
  result.m[3][0] = (left + right) / (left - right);
  result.m[3][1] = (top + bottom) / (bottom - top);
  result.m[3][2] = -nearClip / (farClip - nearClip);
  result.m[3][3] = 1.0f;

  return result;
}

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height,
                             float minDepth, float maxDepth) {
  Matrix4x4 result = {};

  float halfWidth = width * 0.5f;
  float halfHeight = height * 0.5f;
  float depthRange = maxDepth - minDepth;

  result.m[0][0] = halfWidth;
  result.m[1][1] = -halfHeight; // Y 軸を反転（DirectX は左上原点）
  result.m[2][2] = depthRange;
  result.m[3][0] = left + halfWidth;
  result.m[3][1] = top + halfHeight;
  result.m[3][2] = minDepth;
  result.m[3][3] = 1.0f;

  return result;
}

Matrix4x4 Transpose(const Matrix4x4 &m) {
  Matrix4x4 result;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      result.m[i][j] = m.m[j][i];
    }
  }
  return result;
}

Matrix4x4 MakeIdentity4x4() {
  Matrix4x4 result;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        result.m[i][j] = 1.0f;
      } else {
        result.m[i][j] = 0.0f;
      }
    }
  }
  return result;
}

Vector3 Cross(const Vector3 &v1, const Vector3 &v2) {
  Vector3 result;
  result.x = v1.y * v2.z - v1.z * v2.y;
  result.y = v1.z * v2.x - v1.x * v2.z;
  result.z = v1.x * v2.y - v1.y * v2.x;
  return result;
}

Vector3 Project(const Vector3 &v1, const Vector3 &v2) {
  float dot = v1.Dot(v2);
  float lenSq = v2.LengthSquared();
  if (lenSq == 0.0f)
    return {0, 0, 0}; // v2がゼロベクトルなら射影はゼロベクトル
  float scale = dot / lenSq;
  return v2 * scale;
};

Vector3 ClosestPoint(const Vector3 &point, const Segment &segment) {
  Vector3 ab = segment.diff - segment.origin;
  Vector3 ap = point - segment.origin;

  float abLenSq = ab.LengthSquared();
  if (abLenSq == 0.0f)
    return segment.origin; // degenerate segment

  float t = ap.Dot(ab) / abLenSq;

  if (t <= 0.0f)
    return segment.origin;
  if (t >= 1.0f)
    return segment.diff;
  return segment.origin + (segment.diff * t);
};

float Length(const Vector3 &v) {
  return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// bool isCollision(Sphere &s1, Sphere &s2) {
//   float distance = Length(s2.center - s1.center);
//   if (distance <= s1.radius + s2.radius) {
//     s1.color = RED;
//     return true;
//   }
//   return false;
// }

bool IsCollision(const Sphere &sphere, const Plane &plane) {
  float distance = Dot(plane.normal, sphere.center) - plane.distance;

  // 距離の絶対値が半径以下なら衝突している
  return std::abs(distance) <= sphere.radius;
}

Vector3 Perpendicular(const Vector3 &vector) {
  if (vector.x != 0.0f || vector.y != 0.0f) {
    return {-vector.y, vector.x, 0.0f};
  }
  return {0.0f, -vector.z, vector.y};
}

static const int kRowheight = 20;

static const int kColumnWidth = 60;

void Matrix4x4ScreenPrintf(int x, int y, const Matrix4x4 &matrix,
                           const char *lavel) {

  Novice::ScreenPrintf(x, y, lavel);
  for (int row = 0; row < 4; row++) {
    for (int column = 0; column < 4; column++) {
      Novice::ScreenPrintf(x + column * kColumnWidth,
                           y + (row + 1) * kRowheight, "%6.02f",
                           matrix.m[row][column]);
    }
  }
}

void Vector3ScreenPrintf(int x, int y, const Vector3 &vector,
                         const char *label) {
  Novice::ScreenPrintf(x, y, "%.02f", vector.x);
  Novice::ScreenPrintf(x + kColumnWidth, y, "%.02f", vector.y);
  Novice::ScreenPrintf(x + kColumnWidth * 2, y, "%.02f", vector.z);
  Novice::ScreenPrintf(x + kColumnWidth * 3, y, label);
}

void DrawGrid(const Matrix4x4 &viewProjectionMatrix,
              const Matrix4x4 &viewportMatrix) {
  const float kGridHalfWidth = 2.0f;
  const uint32_t kSubdivision = 10;
  const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);
  // 奥から手前への線を順に引いていく（X軸方向）
  for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
    float x = -kGridHalfWidth + xIndex * kGridEvery;

    // ワールド座標系で線の始点と終点を計算
    Vector3 worldStart = {x, 0.0f, -kGridHalfWidth};
    Vector3 worldEnd = {x, 0.0f, kGridHalfWidth};

    // スクリーン座標に変換
    Vector3 clipStart = Transform(worldStart, viewProjectionMatrix);
    Vector3 screenStart = Transform(clipStart, viewportMatrix);

    Vector3 clipEnd = Transform(worldEnd, viewProjectionMatrix);
    Vector3 screenEnd = Transform(clipEnd, viewportMatrix);

    // 線の色（薄いグレー: 0xAAAAAAFF）
    uint32_t color = 0xAAAAAAFF;

    // 原点（x = 0）のときだけ赤色にする
    if (abs(x) < 0.0001f) {
      color = 0xFF0000FF;
    }

    // 線を描画
    Novice::DrawLine(int(screenStart.x), int(screenStart.y), int(screenEnd.x),
                     int(screenEnd.y), color);
  }

  // 左から右への線を順に引いていく（Z軸方向）
  for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
    float z = -kGridHalfWidth + zIndex * kGridEvery;

    Vector3 worldStart = {-kGridHalfWidth, 0.0f, z};
    Vector3 worldEnd = {kGridHalfWidth, 0.0f, z};

    Vector3 clipStart = Transform(worldStart, viewProjectionMatrix);
    Vector3 screenStart = Transform(clipStart, viewportMatrix);

    Vector3 clipEnd = Transform(worldEnd, viewProjectionMatrix);
    Vector3 screenEnd = Transform(clipEnd, viewportMatrix);

    uint32_t color = 0xAAAAAAFF;
    if (abs(z) < 0.0001f) {
      color = 0x0000FFFF; // 原点に交差するZ軸線だけ青色などにしたければ
    }

    Novice::DrawLine(int(screenStart.x), int(screenStart.y), int(screenEnd.x),
                     int(screenEnd.y), color);
  }
}

void DrawSphere(const Sphere &sphere, const Matrix4x4 &viewProjectionMatrix,
                const Matrix4x4 &viewportMatrix, uint32_t color) {
  const uint32_t kSubdivision = 20; // 分割数
  const float kLonEvery =
      2.0f * float(M_PI) / kSubdivision;              // 経度分割1つ分の角度
  const float kLatEvery = float(M_PI) / kSubdivision; // 緯度分割1つ分の角度

  // 緯度の角度を-π/2 〜 π/2にしてループ
  for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
    float lat = -float(M_PI) / 2.0f + kLatEvery * latIndex; // 現在の緯度

    for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
      float lon = lonIndex * kLonEvery; // 現在の経度

      // 緯度経度から位置ベクトルを計算（球面座標系 -> デカルト座標系）
      Vector3 a = {sphere.center.x + sphere.radius * cosf(lat) * cosf(lon),
                   sphere.center.y + sphere.radius * sinf(lat),
                   sphere.center.z + sphere.radius * cosf(lat) * sinf(lon)};

      Vector3 b = {
          sphere.center.x + sphere.radius * cosf(lat + kLatEvery) * cosf(lon),
          sphere.center.y + sphere.radius * sinf(lat + kLatEvery),
          sphere.center.z + sphere.radius * cosf(lat + kLatEvery) * sinf(lon)};

      Vector3 c = {
          sphere.center.x + sphere.radius * cosf(lat) * cosf(lon + kLonEvery),
          sphere.center.y + sphere.radius * sinf(lat),
          sphere.center.z + sphere.radius * cosf(lat) * sinf(lon + kLonEvery)};

      // a, b, cをscreen座標に変換
      Vector3 aScreen =
          Transform(Transform(a, viewProjectionMatrix), viewportMatrix);
      Vector3 bScreen =
          Transform(Transform(b, viewProjectionMatrix), viewportMatrix);
      Vector3 cScreen =
          Transform(Transform(c, viewProjectionMatrix), viewportMatrix);

      // ab, bcで線を描画
      Novice::DrawLine(int(aScreen.x), int(aScreen.y), int(bScreen.x),
                       int(bScreen.y), color);
      Novice::DrawLine(int(aScreen.x), int(aScreen.y), int(cScreen.x),
                       int(cScreen.y), color);
    }
  }
}

void DrawPlane(const Plane &plane, const Matrix4x4 &viewProjectionMatrix,
               const Matrix4x4 &viewportMatrix, uint32_t color) {
  // 法線を正規化（念のため）
  Vector3 normal = Normalize(plane.normal);

  // 平面の中心点（法線方向に distance 分オフセット）
  Vector3 center = Multiply(plane.distance, normal);

  // 平面上の2つの直交ベクトルを生成
  Vector3 right = Normalize(Perpendicular(normal));
  Vector3 forward = Normalize(Cross(normal, right));

  // 平面のサイズ（適宜変更可能）
  float halfSize = 2.0f;

  // 4つの角の座標（時計回りに）
  Vector3 corners[4];
  corners[0] = Add(center, Add(Multiply(+halfSize, right),
                               Multiply(+halfSize, forward))); // 右前
  corners[1] = Add(center, Add(Multiply(-halfSize, right),
                               Multiply(+halfSize, forward))); // 左前
  corners[2] = Add(center, Add(Multiply(-halfSize, right),
                               Multiply(-halfSize, forward))); // 左後
  corners[3] = Add(center, Add(Multiply(+halfSize, right),
                               Multiply(-halfSize, forward))); // 右後

  // 各頂点をスクリーン座標へ変換
  for (int i = 0; i < 4; ++i) {
    corners[i] =
        Transform(Transform(corners[i], viewProjectionMatrix), viewportMatrix);
  }

  // 線を描画（4辺）
  for (int i = 0; i < 4; ++i) {
    const Vector3 &p1 = corners[i];
    const Vector3 &p2 = corners[(i + 1) % 4];

    int x1 = static_cast<int>(p1.x);
    int y1 = static_cast<int>(p1.y);
    int x2 = static_cast<int>(p2.x);
    int y2 = static_cast<int>(p2.y);

    Novice::DrawLine(x1, y1, x2, y2, color);
  }
}

void DrawSegment(const Segment &segment, const Matrix4x4 &viewProjectionMatrix,
                 const Matrix4x4 &viewportMatrix, uint32_t color) {
  // 始点と終点をスクリーン座標に変換
  Vector3 startScreen = Transform(
      Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
  Vector3 endScreen =
      Transform(Transform(segment.diff, viewProjectionMatrix), viewportMatrix);
  // 線を描画
  Novice::DrawLine(int(startScreen.x), int(startScreen.y), int(endScreen.x),
                   int(endScreen.y), color);
}

bool IsCollision(const Segment &line, const Plane &plane) {
  // まず垂直判定を行うために、法線と線の方向の内積を求める
  float dot = Dot(plane.normal, line.diff);

  // 垂直＝平行であるので、衝突しているはずがない
  if (dot == 0.0f) {
    return false;
  }

  // t を求める
  float t = (plane.distance - Dot(line.origin, plane.normal)) / dot;

  // t の値と線分の範囲によって衝突しているかを判断する
  return t >= 0.0f && t <= 1.0f;
}

#pragma endregion

const int kWindowWidth = 1280;
const int kWindowHeight = 720;

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, 1280, 720);

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

#pragma region カメラ用変数
  const float rotateSensitivity = 0.003f;
  const float zoomSensitivity = 0.05f;
  const float panSensitivity = 0.01f;
  static int prevMouseX = 0, prevMouseY = 0;
  int mouseX, mouseY;

  // カメラ用変数（グローバル or クラス内）
  Vector3 target = {0.0f, 0.0f, 0.0f};
  float distance = 10.0f;
  Vector3 cameraTranslate{0.0f, 1.9f, -6.49f};
  Vector3 cameraRotate{0.26f, 0.0f, 0.0f};

  Camera camera = {cameraTranslate, 0.0f, 0.0f};

#pragma endregion

  Plane plane = {Normalize({0.0f, 1.0f, 1.0f}), 0.5f}; // 斜め上向き

  Segment segment = {
      {0.0f, 0.0f, 0.0f},   // 始点
      {1.0f, 1.0f, 1.0f} // 終点
  };

  unsigned int color = WHITE; // 球の色

  // ウィンドウの×ボタンが押されるまでループ
  while (Novice::ProcessMessage() == 0) {
    // フレームの開始
    Novice::BeginFrame();

    // キー入力を受け取る
    memcpy(preKeys, keys, 256);
    Novice::GetHitKeyStateAll(keys);

    ///
    /// ↓更新処理ここから
    ///

#pragma region カメラ操作+renderingPipeline
    // 入力
    Novice::GetMousePosition(&mouseX, &mouseY);

    // 回転（右クリックドラッグ）
    if (!ImGui::GetIO().WantCaptureMouse && Novice::IsPressMouse(1)) {
      int deltaX = mouseX - prevMouseX;
      int deltaY = mouseY - prevMouseY;

      cameraRotate.y += deltaX * rotateSensitivity;
      cameraRotate.x -= deltaY * rotateSensitivity;

      // 上下の回転制限（Blenderっぽく ±85°）
      cameraRotate.x = std::clamp(cameraRotate.x, -1.48f, 1.48f);
    }

    // ズーム（マウスホイール）
    int wheel = Novice::GetWheel();
    if (wheel != 0) {
      distance -= wheel * zoomSensitivity;
      distance = std::max(distance, 1.0f); // 最小距離
    }

    // パン（ホイールクリックドラッグ）
    if (!ImGui::GetIO().WantCaptureMouse && Novice::IsPressMouse(2)) {
      int deltaX = mouseX - prevMouseX;
      int deltaY = mouseY - prevMouseY;

      // カメラの右ベクトルと上ベクトルでターゲットを移動
      Vector3 right = {cosf(cameraRotate.y), 0, -sinf(cameraRotate.y)};
      Vector3 up = {0, 1, 0};

      target.x -= (right.x * deltaX + up.x * -deltaY) * panSensitivity;
      target.y -= (right.y * deltaX + up.y * -deltaY) * panSensitivity;
      target.z -= (right.z * deltaX + up.z * -deltaY) * panSensitivity;
    }

    prevMouseX = mouseX;
    prevMouseY = mouseY;

    // カメラ座標計算
    cameraTranslate.x =
        target.x + distance * cosf(cameraRotate.x) * sinf(cameraRotate.y);
    cameraTranslate.y = target.y + distance * sinf(cameraRotate.x);
    cameraTranslate.z =
        target.z + distance * cosf(cameraRotate.x) * cosf(cameraRotate.y);

    // 行列計算
    Matrix4x4 cameraMatrix =
        MakeAffineMatrix({1.0f, 1.0f, 1.0f}, cameraRotate, cameraTranslate);
    Matrix4x4 viewMatrix = Inverse(cameraMatrix);

    Matrix4x4 worldMatrix = MakeIdentity4x4();
    Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(
        0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
    Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
    Matrix4x4 viewportMatrix = MakeViewportMatrix(
        0.0f, 0.0f, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);
#pragma endregion

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    DrawGrid(viewProjectionMatrix, viewportMatrix);
    DrawPlane(plane, viewProjectionMatrix, viewportMatrix, color);
    if (IsCollision(segment, plane)) {
      DrawSegment(segment, viewProjectionMatrix, viewportMatrix, RED);
    } else {
      DrawSegment(segment, viewProjectionMatrix, viewportMatrix, WHITE);
    }

    ImGui::Begin("Window");

    // Plane.Normal
    ImGui::DragFloat3("Plane.Normal", &plane.normal.x, 0.01f);

    // Plane.Distance
    ImGui::DragFloat("Plane.Distance", &plane.distance, 0.01f);

    ImGui::DragFloat3("Origin", &segment.origin.x, 0.1f);

    ImGui::DragFloat3("Diff", &segment.diff.x, 0.1f);

    ImGui::End();

    ///
    /// ↑描画処理ここまで
    ///

    // フレームの終了
    Novice::EndFrame();

    // ESCキーが押されたらループを抜ける
    if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
      break;
    }
  }

  // ライブラリの終了
  Novice::Finalize();
  return 0;
}
