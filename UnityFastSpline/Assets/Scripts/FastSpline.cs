using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// TODO:
// - time calc should probably weight two adjacent lengths
// - transform handle 'soft' mode, pull nearby verts with falloff
// - nearest point
// - distance walk

#if UNITY_EDITOR
using UnityEditor;

[CustomEditor(typeof(FastSpline))]
public class FastSplineEditor
    : Editor
{
    private int selectedIndex = -1;

    private bool preview;
    private bool previewAuto;
    private float previewTime;

    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        var self = target as FastSpline;

        GUILayout.Label("");

        GUILayout.Label(string.Format("Points: {0}", self.points.Count));

        GUILayout.Label(string.Format("Approximate Arc Length: {0}", self.CalculateApproximateArcLength()));
        GUILayout.Label(string.Format("Longest segment: {0}", self.CalculateLongestSegment()));

        var lengthTotal = self.CalculateApproximateArcLength();
        var lengthLongest = self.CalculateLongestSegment();
        var longestRatio = lengthLongest / lengthTotal;

        if (longestRatio > 0.15f)
            EditorGUILayout.HelpBox("single point quite far outlier", MessageType.Warning);

        for (int i = 0; i < self.lengths.Count; ++i)
            GUILayout.Label(string.Format("Length: {0}-{1}: {2}", i, i + 1, self.lengths[i]));
        for (int i = 0; i < self.times.Count; ++i)
            GUILayout.Label(string.Format("Times: {0}: {1}", i, self.times[i]));

        GUILayout.Label("");

        preview = GUILayout.Toggle(preview, "Preview");
        previewTime = EditorGUILayout.Slider("Preview Time", previewTime, 0.0f, 1.0f);

        if (GUILayout.Button("TEST PB SPLINE"))
            PaulBourkeSplineTest.TestSpline();
    }

    public void OnSceneGUI()
    {
        var self = target as FastSpline;

        for (int i = 0; i < self.points.Count; ++i)
        {
            var p0 = self.points[i];
            var w0 = self.transform.TransformPoint(p0);

            if (selectedIndex == i)
            {
                //Handles.color = Color.blue;
                //Handles.Button(w0, Quaternion.identity, 0.25f, 0.25f, Handles.RectangleHandleCap);

                var refPosition = p0;
                var refRotation = Quaternion.identity;

                refPosition = self.transform.TransformPoint(refPosition);

                Handles.TransformHandle(ref refPosition, ref refRotation);

                refPosition = self.transform.InverseTransformPoint(refPosition);

                var dirty = self.points[i] != refPosition;
                if (dirty)
                {
                    Undo.RecordObject(self, "Move spline point");
                    self.points[i] = refPosition;
                    self.RecalculateInternals();
                }
            }
            else
            {
                Handles.color = Color.white;

                if (Handles.Button(w0, Quaternion.identity, 0.25f, 0.25f, Handles.RectangleHandleCap))
                    selectedIndex = i;
            }
        }

        if (preview)
        {
            var p0 = self.CalculatePosition(previewTime);
            var w0 = self.transform.TransformPoint(p0);

            Handles.color = Color.red;
            Handles.DrawSphere(-1, w0, Quaternion.identity, 0.5f);

            var pn0 = self.CalculatePositionNormalized(previewTime);
            var wn0 = self.transform.TransformPoint(pn0);

            Handles.color = Color.green;
            Handles.DrawCube(-1, wn0, Quaternion.identity, 0.5f);
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


// points: N points describing spline
// lengths: linear lengths between each segment

public class FastSpline
    : MonoBehaviour
{
    public List<Vector3> points;

    [System.NonSerialized]
    public List<float> times;

    [System.NonSerialized]
    public List<float> lengths;

    [System.NonSerialized]
    public float lengthsTotal;

    public FastSpline()
    {
        if (points == null)
        {
            points = new List<Vector3>();
            points.Add(Vector3.zero);
            points.Add(Vector3.zero);
            points.Add(Vector3.zero);

            RecalculateInternals();
        }
    }

    public void OnValidate()
    {
        RecalculateInternals();
    }

    public void RecalculateInternals()
    {
        lengths = new List<float>(points.Count - 1);
        lengthsTotal = 0.0f;

        for (int i = 0; i < points.Count - 1; ++i)
        {
            var p0 = points[i];
            var p1 = points[i + 1];
            var ofs = p1 - p0;
            var len = ofs.sqrMagnitude;
            if (len > 0.00001f)
                len = Mathf.Sqrt(len);

            lengths.Add(len);
            lengthsTotal += len;
        }

        times = new List<float>(points.Count);

        var step = 1.0f / (points.Count - 1);
        var time = 0.0f;

        times.Add(0.0f);

        for (int i = 0; i < lengths.Count; ++i)
        {
            var expected = step;
            var actual = lengths[i];
            var ratio = expected / actual;

            time += step * ratio * lengthsTotal;

            times.Add(time);
        }

        var timeNormalize = 1.0f / times[times.Count - 1];
        for (int i = 0; i < times.Count; ++i)
            times[i] = times[i] * timeNormalize;
    }

    public int CalculateIndex(float t)
    {
        var at = t * lengthsTotal;

        for (int i = 0; i < lengths.Count; ++i)
        {
            var length = lengths[i];

            at -= length;

            if (at <= 0.0f)
                return i;
        }

        return points.Count - 1;
    }

    public Vector3 CalculatePosition(float t)
    {
        var local = t * (points.Count - 1);
        var index = (int)local;
        var time = local - (float)index;

        var a = index - 1;
        var b = index;
        var c = index + 1;
        var d = index + 2;

        a = Mathf.Max(a, 0);
        b = Mathf.Min(b, points.Count - 1);
        c = Mathf.Min(c, points.Count - 1);
        d = Mathf.Min(d, points.Count - 1);

        var p0 = points[a];
        var p1 = points[b];
        var p2 = points[c];
        var p3 = points[d];

        var t1 = time;
        var t2 = time * time;
        var t3 = time * time * time;

        var result = 0.5f * ((2.0f * p1) + (-p0 + p2) * t1 + (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 + (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);

        return result;
    }

    public float CalculateTimeNormalized(float t)
    {
        var local = t * (times.Count - 1);
        var index = (int)local;
        var time = local - (float)index;

        var a = index - 1;
        var b = index;
        var c = index + 1;
        var d = index + 2;

        a = Mathf.Max(a, 0);
        b = Mathf.Min(b, times.Count - 1);
        c = Mathf.Min(c, times.Count - 1);
        d = Mathf.Min(d, times.Count - 1);

        var p0 = times[a];
        var p1 = times[b];
        var p2 = times[c];
        var p3 = times[d];

        var t1 = time;
        var t2 = time * time;
        var t3 = time * time * time;

        var result = 0.5f * ((2.0f * p1) + (-p0 + p2) * t1 + (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 + (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
        var rescaled = result;

        if (float.IsNaN(rescaled))
            rescaled = 0.0f;
        rescaled = Mathf.Clamp01(rescaled);

        return rescaled;
    }

    public Vector3 CalculatePositionNormalized(float t)
    {
        var resultTime = CalculateTimeNormalized(t);
        var resultPosition = CalculatePosition(resultTime);

        return resultPosition;
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
            var lsq = ofs.sqrMagnitude;
            if (lsq > Mathf.Epsilon)
                lsq = Mathf.Sqrt(lsq);

            sum += lsq;
        }

        return sum;
    }

    public float CalculateLongestSegment()
    {
        var longest = 0.0f;

        for (int i = 1; i < lengths.Count; ++i)
        {
            var length = lengths[i];
            if (length > longest)
                longest = length;
        }

        return longest;
    }

    public Vector3 Project(Vector3 p)
    {
        var iterations = 8;
        var time = 0.0f;

        return Vector3.zero;
    }

#if UNITY_EDITOR
    public void OnDrawGizmos()
    {
        var steps = 100;

        for (int i = 1; i < steps; ++i)
        {
            var t0 = 1.0f / (float)steps * (float)(i - 1);
            var t1 = 1.0f / (float)steps * (float)(i);
            var p0 = CalculatePosition(t0);
            var p1 = CalculatePosition(t1);
            var w0 = transform.TransformPoint(p0);
            var w1 = transform.TransformPoint(p1);

            Gizmos.color = Color.Lerp(Color.white, Color.clear, 0.5f);
            Gizmos.DrawLine(w0, w1);
        }
    }

    public void OnDrawGizmosSelected()
    {
        var steps = 100;

        for (int i = 1; i < steps; ++i)
        {
            var t0 = 1.0f / (float)steps * (float)(i - 1);
            var t1 = 1.0f / (float)steps * (float)(i);
            var p0 = CalculatePosition(t0);
            var p1 = CalculatePosition(t1);
            var w0 = transform.TransformPoint(p0);
            var w1 = transform.TransformPoint(p1);

            Gizmos.color = Color.Lerp(Color.white, Color.clear, 0.5f);
            Gizmos.DrawLine(w0, w1);
        }

        for (int i = 1; i < points.Count; ++i)
        {
            var p0 = points[i - 1];
            var p1 = points[i];
            var w0 = transform.TransformPoint(p0);
            var w1 = transform.TransformPoint(p1);

            Gizmos.color = Color.Lerp(Color.white, Color.blue, 0.5f);
            Gizmos.DrawLine(w0, w1);
        }

        for (int i = 0; i < points.Count; ++i)
        {
            var t = 1.0f / (points.Count - 1) * i;
            var p = CalculatePositionNormalized(t);
            var w = transform.TransformPoint(p);

            Gizmos.color = Color.Lerp(Color.red, Color.blue, t);
            Gizmos.DrawSphere(w, 0.25f);
        }

        for (int i = 0; i < times.Count; ++i)
        {
            var a = points[0];
            var b = points[points.Count - 1];
            var t = times[i];
            var p = a + (b - a) * t;
            var w = transform.TransformPoint(p);

            Gizmos.color = Color.Lerp(Color.green, Color.white, t);
            Gizmos.DrawCube(w, Vector3.one * 0.25f);
        }

        for (int i = 0; i < 50; ++i)
        {
            var a = points[0];
            var b = points[points.Count - 1];
            var t = CalculateTimeNormalized(1.0f / 50.0f * i);
            var p = a + (b - a) * t;
            var w = transform.TransformPoint(p);

            Gizmos.color = Color.Lerp(Color.blue, Color.blue, t);
            Gizmos.DrawCube(w-Vector3.up * 0.5f, Vector3.one * 0.10f);
        }
    }
#endif
}

// ref: http://paulbourke.net/geometry/spline/

public class PaulBourkeSplineTest
{
    // This returns the point "output" on the spline curve.
    // The parameter "v" indicates the position, it ranges from 0 to n-t+2
    
    public static Vector3 SplinePoint(int[] u, int n, int t, float v, Vector3[] control)
    {
        Vector3 output = Vector3.zero;

        int k;
        float b;

        for (k = 0; k <= n; k++)
        {
            b = SplineBlend(k, t, u, v);
            output.x += control[k].x * b;
            output.y += control[k].y * b;
            output.z += control[k].z * b;
        }

        return output;
    }

    // Calculate the blending value, this is done recursively.
    // If the numerator and denominator are 0 the expression is 0.
    // If the deonimator is 0 the expression is 0

    public static float SplineBlend(int k, int t, int[] u, float v)
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

    public static void SplineKnots(int[] u, int n, int t)
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

    public static void SplineCurve(Vector3[] inp, int n, int[] knots, int t, Vector3[] outp, int res)
    {
        int i;
        float interval, increment;

        interval = 0;
        increment = (n - t + 2) / (float)(res - 1);
        for (i = 0; i < res - 1; i++)
        {
            outp[i] = SplinePoint(knots, n, t, interval, inp);
            interval += increment;
        }
        outp[res - 1] = inp[n];
    }


    // Example of how to call the spline functions
    // Basically one needs to create the control points, then compute
    // the knot positions, then calculate points along the curve.

    public static int N = 3;
    public static int T = 3;
    public static int RESOLUTION = 200;
    public static int[] knots = new int[N + T + 1];
    public static Vector3[] outp = new Vector3[RESOLUTION];

    public static void TestSpline()
    {
        var inp = new Vector3[]
        {
            new Vector3(0.0f, 0.0f, 0.0f),
            new Vector3(1.0f, 0.0f, 3.0f),
            new Vector3(2.0f, 0.0f, 1.0f),
            new Vector3(4.0f, 0.0f, 4.0f),
        };

        SplineKnots(knots, N, T);
        SplineCurve(inp, N, knots, T, outp, RESOLUTION);

        for (int i = 0; i < RESOLUTION; ++i)
            Debug.LogFormat("{0}: {1},{2},{3}", i, outp[i].x, outp[i].y, outp[i].z);
    }
}
