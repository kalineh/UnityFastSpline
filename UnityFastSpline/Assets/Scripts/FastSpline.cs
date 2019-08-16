using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/*
// ref: http://paulbourke.net/geometry/spline/

// This returns the point "output" on the spline curve.
// The parameter "v" indicates the position, it ranges from 0 to n-t+2

void SplinePoint(int* u, int n, int t, float v, XYZ* control, XYZ* output)
{
    int k;
    float b;

    output->x = 0;
    output->y = 0;
    output->z = 0;

    for (k = 0; k <= n; k++)
    {
        b = SplineBlend(k, t, u, v);
        output->x += control[k].x * b;
        output->y += control[k].y * b;
        output->z += control[k].z * b;
    }
}

// Calculate the blending value, this is done recursively.
// If the numerator and denominator are 0 the expression is 0.
// If the deonimator is 0 the expression is 0

float SplineBlend(int k, int t, int* u, float v)
{
    float value;

    if (t == 1)
    {
        if ((u[k] <= v) && (v < u[k + 1]))
            value = 1;
        else
            value = 0;
    }
    else
    {
        if ((u[k + t - 1] == u[k]) && (u[k + t] == u[k + 1]))
            value = 0;
        else if (u[k + t - 1] == u[k])
            value = (u[k + t] - v) / (u[k + t] - u[k + 1]) * SplineBlend(k + 1, t - 1, u, v);
        else if (u[k + t] == u[k + 1])
            value = (v - u[k]) / (u[k + t - 1] - u[k]) * SplineBlend(k, t - 1, u, v);
        else
            value = (v - u[k]) / (u[k + t - 1] - u[k]) * SplineBlend(k, t - 1, u, v) +
                    (u[k + t] - v) / (u[k + t] - u[k + 1]) * SplineBlend(k + 1, t - 1, u, v);
    }
    return (value);
}


// The positions of the subintervals of v and breakpoints, the position
// on the curve are called knots. Breakpoints can be uniformly defined
// by setting u[j] = j, a more useful series of breakpoints are defined
// by the function below. This set of breakpoints localises changes to
// the vicinity of the control point being modified.

void SplineKnots(int* u, int n, int t)
{
    int j;

    for (j = 0; j <= n + t; j++)
    {
        if (j < t)
            u[j] = 0;
        else if (j <= n)
            u[j] = j - t + 1;
        else if (j > n)
            u[j] = n - t + 2;
    }
}

// Create all the points along a spline curve
// Control points "inp", "n" of them.
// Knots "knots", degree "t".
// Ouput curve "outp", "res" of them.

void SplineCurve(XYZ* inp, int n, int* knots, int t, XYZ* outp, int res)
{
    int i;
    float interval, increment;

    interval = 0;
    increment = (n - t + 2) / (float)(res - 1);
    for (i = 0; i < res - 1; i++)
    {
        SplinePoint(knots, n, t, interval, inp, &(outp[i]));
        interval += increment;
    }
    outp[res - 1] = inp[n];
}


// Example of how to call the spline functions
// Basically one needs to create the control points, then compute
// the knot positions, then calculate points along the curve.

#define N 3
XYZ inp[N + 1] = { 0.0, 0.0, 0.0, 1.0, 0.0, 3.0, 2.0, 0.0, 1.0, 4.0, 0.0, 4.0 };
#define T 3
int knots[N + T + 1];
#define RESOLUTION 200
XYZ outp[RESOLUTION];

int main(int argc, char** argv)
{
    int i;

    SplineKnots(knots, N, T);
    SplineCurve(inp, N, knots, T, outp, RESOLUTION);

    // Display the curve, in this case in OOGL format for GeomView
    printf("LIST\n");
    printf("{ = SKEL\n");
    printf("%d %d\n", RESOLUTION, RESOLUTION - 1);
    for (i = 0; i < RESOLUTION; i++)
        printf("%g %g %g\n", outp[i].x, outp[i].y, outp[i].z);
    for (i = 0; i < RESOLUTION - 1; i++)
        printf("2 %d %d 1 1 1 1\n", i, i + 1);
    printf("}\n");

    // The axes
    printf("{ = SKEL 3 2  0 0 4  0 0 0  4 0 0  2 0 1 0 0 1 1 2 1 2 0 0 1 1 }\n");

    // Control point polygon
    printf("{ = SKEL\n");
    printf("%d %d\n", N + 1, N);
    for (i = 0; i <= N; i++)
        printf("%g %g %g\n", inp[i].x, inp[i].y, inp[i].z);
    for (i = 0; i < N; i++)
        printf("2 %d %d 0 1 0 1\n", i, i + 1);
    printf("}\n");

}
*/

#if UNITY_EDITOR
using UnityEditor;

[CustomEditor(typeof(FastSpline))]
public class FastSplineEditor
    : Editor
{
    private int selectedIndex = -1;

    private bool preview;
    private float previewTime;

    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        var self = target as FastSpline;

        GUILayout.Label("");

        GUILayout.Label(string.Format("Points: {0}", self.points.Count));
        GUILayout.Label(string.Format("Approximate Arc Length: {0}", self.CalculateApproximateArcLength()));

        GUILayout.Label("");

        preview = GUILayout.Toggle(preview, "Preview");
        previewTime = EditorGUILayout.Slider("Preview Time", previewTime, 0.0f, 1.0f);
    }

    public void OnSceneGUI()
    {
        var self = target as FastSpline;
        var steps = 10;

        for (int i = 1; i < steps; ++i)
        {
            var t0 = 1.0f / (float)steps * (float)(i - 1);
            var t1 = 1.0f / (float)steps * (float)(i);
            var p0 = self.CalculatePosition(t0);
            var p1 = self.CalculatePosition(t1);
            var w0 = self.transform.TransformPoint(p0);
            var w1 = self.transform.TransformPoint(p1);

            Handles.DrawLine(w0, w1);
        }

        for (int i = 0; i < self.points.Count; ++i)
        {
            var p0 = self.points[i];
            var w0 = self.transform.TransformPoint(p0);

            if (selectedIndex == i)
            {
                Handles.color = Color.blue;
                Handles.Button(w0, Quaternion.identity, 0.5f, 0.5f, Handles.RectangleHandleCap);

                var refPosition = p0;
                var refRotation = Quaternion.identity;

                refPosition = self.transform.TransformPoint(refPosition);

                Handles.TransformHandle(ref refPosition, ref refRotation);

                refPosition = self.transform.InverseTransformPoint(refPosition);

                self.points[i] = refPosition;
            }
            else
            {
                Handles.color = Color.white;

                if (Handles.Button(w0, Quaternion.identity, 0.5f, 0.5f, Handles.RectangleHandleCap))
                    selectedIndex = i;
            }
        }

        if (preview)
        {
            var p0 = self.CalculatePosition(previewTime);
            var w0 = self.transform.TransformPoint(p0);

            Handles.color = Color.red;
            Handles.DrawSphere(-1, w0, Quaternion.identity, 0.5f);
        }
    }
}
#endif


public static class FastSplineMath
{
    public static Vector3 QuadraticInterpolate(List<Vector3> points, int index, float t)
    {
        var count = points.Count;

        var index0 = Mathf.Clamp(index + 0, 0, count - 1);
        var index1 = Mathf.Clamp(index + 1, 0, count - 1);
        var index2 = Mathf.Clamp(index + 2, 0, count - 1);

        var p0 = points[index0];
        var p1 = points[index1];
        var p2 = points[index2];

        return QuadraticInterpolate(p0, p1, p2, t);
    }

    public static Vector3 QuadraticInterpolate(Vector3 p0, Vector3 p1, Vector3 p2, float t)
    {
        var tt = t * t;
        var u = 1.0f - t;
        var uu = u * u;
        var p = p0 * uu + p1 * 2.0f * u * t + p2 * tt;

        return p;
    }

    public static Vector3 CubicInterpolate(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, float t)
    {
        var u = 1.0f - t;
        var tt = t * t;
        var uu = u * u;
        var uuu = uu * u;
        var ttt = tt * t;
        var p = uuu * p0;

        p += 3.0f * uu * t * p1;
        p += 3.0f* u * tt * p2;
        p += ttt * p3;

        return p;
    }
}

public class FastSpline
    : MonoBehaviour
{
    public List<Vector3> points;

    private List<Vector3> knots;

    public FastSpline()
    {
        if (points == null)
        {
            points = new List<Vector3>();
            points.Add(Vector3.zero);
            points.Add(Vector3.zero);
            points.Add(Vector3.zero);
        }
    }

    public Vector3 CalculatePosition(float t)
    {
        var index = (int)(t * (float)points.Count);
        var time = (t / (float)points.Count) * 3.0f;
        var point = FastSplineMath.QuadraticInterpolate(points, index, t);

        return point;
    }

    public float CalculateApproximateArcLength()
    {
        var sum = 0.0f;
        var segments = 100;

        for (int i = 1; i < segments; ++i)
        {
            var t0 = 1.0f / (float)segments * (float)(i - 1);
            var t1 = 1.0f / (float)segments * (float)(i);
            var p0 = CalculatePosition(t0);
            var p1 = CalculatePosition(t1);
            var ofs = (p1 - p0);

            sum += ofs.sqrMagnitude; 
        }

        if (sum > 0.00001f)
            sum = Mathf.Sqrt(sum);

        return sum;
    }

#if UNITY_EDITOR
    public void OnDrawGizmos()
    {
        var steps = 10;

        for (int i = 1; i < steps; ++i)
        {
            var t0 = 1.0f / (float)steps * (float)(i - 1);
            var t1 = 1.0f / (float)steps * (float)(i);
            var p0 = CalculatePosition(t0);
            var p1 = CalculatePosition(t1);
            var w0 = transform.TransformPoint(p0);
            var w1 = transform.TransformPoint(p1);

            Gizmos.color = Color.Lerp(Color.blue, Color.clear, 0.5f);
            Gizmos.DrawLine(w0, w1);
        }
    }
#endif
}
