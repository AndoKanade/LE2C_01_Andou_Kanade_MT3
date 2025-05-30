#include <Novice.h>
#define _USE_MATH_DEFINES
#include <assert.h>
#include <cmath>
#include <imgui.h>
#include <math.h>

const char kWindowTitle[] = "LE2C_01_アンドウ_カナデ_MT3";

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

#pragma region 関数

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
  return segment.origin + ab * t;
};

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

  Vector3 cameraTranslate{0.0f, 1.9f, -6.49f};
  Vector3 cameraRotate{0.26f, 0.0f, 0.0f};

  Segment segment{{-2.0f, -1.0f, 0.0f}, {3.0f, 2.0f, 2.0f}};
  Vector3 point{-1.5f, 0.6f, 0.6f};

  Vector3 project = Project(Subtract(point, segment.origin), segment.diff);

  Vector3 clossestPoint = ClosestPoint(point, segment);

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

    Matrix4x4 cameraMatrix =
        MakeAffineMatrix({1.0f, 1.0f, 1.0f}, cameraRotate, cameraTranslate);
    Matrix4x4 viewMatrix = Inverse(cameraMatrix);

    // ワールド行列は単純な単位行列（またはオブジェクトの変換）
    Matrix4x4 worldMatrix = MakeIdentity4x4();

    // プロジェクションと合成
    Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(
        0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
    Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);

    // ビューポート行列
    Matrix4x4 viewportMatrix = MakeViewportMatrix(
        0.0f, 0.0f, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

    Sphere pointSphere{point, 0.01f};
    Sphere clossestSphire{clossestPoint, 0.01f};

    Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix),
                              viewportMatrix);
    Vector3 end = Transform(
        Transform(Add(segment.origin, segment.diff), viewProjectionMatrix),
        viewportMatrix);

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    DrawSphere(pointSphere, viewProjectionMatrix, viewportMatrix, RED);
    DrawSphere(clossestSphire, viewProjectionMatrix, viewportMatrix, BLACK);
    Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);

    ImGui::InputFloat3("Project", &project.x, "%0.3f",
                       ImGuiInputTextFlags_ReadOnly);

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
