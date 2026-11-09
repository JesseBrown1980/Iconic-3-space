//! transfer_enum: exact enumeration of gauge-fixed lattice configurations
//! for a finite group given by its multiplication table.
//!
//! Integer only. No float anywhere. Input and output are plain text.
//!
//! Input (whitespace-separated integers):
//!   n k np nobs
//!   mult[n*n]   (row-major, index 0 = identity)
//!   inv[n]
//!   ta[n] tb[n] (trace of element g = ta[g] + tb[g]*sqrt(d))
//!   ntr tid[n]  (trace-class id of each element, 0..ntr)
//!   np words, then nobs words; a word is: len, then len pairs (free_link_index, sign)
//!
//! Output:
//!   TOTAL <configs>
//!   ACT <A> <B> <count>            histogram of the action trace sum T = A + B*sqrt(d)
//!   OBS <o> <A> <B> <tid> <count>  joint histogram of (T, trace class of observable o)

use std::env;
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};

type Word = Vec<(usize, bool)>;

struct Input {
    n: usize,
    k: usize,
    mult: Vec<usize>,
    inv: Vec<usize>,
    ta: Vec<i64>,
    tb: Vec<i64>,
    ntr: usize,
    tid: Vec<usize>,
    plaq: Vec<Word>,
    obs: Vec<Word>,
}

struct Reader {
    v: Vec<i64>,
    p: usize,
}

impl Reader {
    fn next(&mut self) -> i64 {
        let x = self.v[self.p];
        self.p += 1;
        x
    }
    fn next_usize(&mut self) -> usize {
        let x = self.next();
        assert!(x >= 0, "negative where an index was expected");
        x as usize
    }
    fn word(&mut self, k: usize) -> Word {
        let len = self.next_usize();
        let mut w = Vec::with_capacity(len);
        for _ in 0..len {
            let idx = self.next_usize();
            let s = self.next();
            assert!(idx < k, "free link index out of range");
            assert!(s == 1 || s == -1, "sign must be +1 or -1");
            w.push((idx, s < 0));
        }
        w
    }
}

fn parse(text: &str) -> Input {
    let v: Vec<i64> = text
        .split_ascii_whitespace()
        .map(|t| t.parse::<i64>().expect("integer token"))
        .collect();
    let mut r = Reader { v, p: 0 };
    let n = r.next_usize();
    let k = r.next_usize();
    let np = r.next_usize();
    let no = r.next_usize();
    let mult: Vec<usize> = (0..n * n).map(|_| r.next_usize()).collect();
    let inv: Vec<usize> = (0..n).map(|_| r.next_usize()).collect();
    let ta: Vec<i64> = (0..n).map(|_| r.next()).collect();
    let tb: Vec<i64> = (0..n).map(|_| r.next()).collect();
    let ntr = r.next_usize();
    let tid: Vec<usize> = (0..n).map(|_| r.next_usize()).collect();
    let plaq: Vec<Word> = (0..np).map(|_| r.word(k)).collect();
    let obs: Vec<Word> = (0..no).map(|_| r.word(k)).collect();
    assert!(r.p == r.v.len(), "trailing tokens in input");
    assert!(mult.iter().all(|&x| x < n), "mult entry out of range");
    assert!(inv.iter().all(|&x| x < n), "inv entry out of range");
    assert!(tid.iter().all(|&x| x < ntr), "tid out of range");
    for g in 0..n {
        assert!(mult[g] == g && mult[g * n] == g, "index 0 is not the identity");
        assert!(mult[g * n + inv[g]] == 0, "inv table wrong");
    }
    Input {
        n,
        k,
        mult,
        inv,
        ta,
        tb,
        ntr,
        tid,
        plaq,
        obs,
    }
}

struct Plan {
    plaq_at: Vec<Vec<usize>>,
    obs_at: Vec<Vec<usize>>,
    const_a: i64,
    const_b: i64,
    const_obs: Vec<usize>,
    width: usize,
    off: i64,
    prefix: usize,
}

fn plan(inp: &Input) -> Plan {
    let mut plaq_at = vec![Vec::new(); inp.k];
    let mut obs_at = vec![Vec::new(); inp.k];
    let mut const_a = 0i64;
    let mut const_b = 0i64;
    let mut const_obs = Vec::new();
    for (p, w) in inp.plaq.iter().enumerate() {
        match w.iter().map(|&(i, _)| i).max() {
            Some(lv) => plaq_at[lv].push(p),
            None => {
                const_a += inp.ta[0];
                const_b += inp.tb[0];
            }
        }
    }
    for (o, w) in inp.obs.iter().enumerate() {
        match w.iter().map(|&(i, _)| i).max() {
            Some(lv) => obs_at[lv].push(o),
            None => const_obs.push(o),
        }
    }
    let np = inp.plaq.len();
    let maxa: i64 = inp.ta.iter().map(|x| x.abs()).max().unwrap_or(0);
    let maxb: i64 = inp.tb.iter().map(|x| x.abs()).max().unwrap_or(0);
    let bound = (maxa.max(maxb)) * (np as i64);
    let width = (2 * bound + 1) as usize;
    let prefix = inp.k.min(2);
    Plan {
        plaq_at,
        obs_at,
        const_a,
        const_b,
        const_obs,
        width,
        off: bound,
        prefix,
    }
}

struct Walker<'a> {
    inp: &'a Input,
    pl: &'a Plan,
    digit: Vec<usize>,
    psa: Vec<i64>,
    psb: Vec<i64>,
    otid: Vec<usize>,
    act: Vec<u64>,
    obsh: Vec<u64>,
}

impl<'a> Walker<'a> {
    fn new(inp: &'a Input, pl: &'a Plan) -> Self {
        let w = pl.width;
        let mut otid = vec![0usize; inp.obs.len()];
        for &o in &pl.const_obs {
            otid[o] = inp.tid[0];
        }
        Walker {
            inp,
            pl,
            digit: vec![0; inp.k],
            psa: vec![0; inp.k],
            psb: vec![0; inp.k],
            otid,
            act: vec![0; w * w],
            obsh: vec![0; inp.obs.len() * w * w * inp.ntr],
        }
    }

    fn hol(&self, word: &[(usize, bool)]) -> usize {
        let n = self.inp.n;
        let mut g = 0usize;
        for &(idx, inverse) in word {
            let mut h = self.digit[idx];
            if inverse {
                h = self.inp.inv[h];
            }
            g = self.inp.mult[g * n + h];
        }
        g
    }

    fn eval(&mut self, j: usize) {
        let (mut a, mut b) = if j == 0 {
            (self.pl.const_a, self.pl.const_b)
        } else {
            (self.psa[j - 1], self.psb[j - 1])
        };
        let pl = self.pl;
        let inp = self.inp;
        for &p in &pl.plaq_at[j] {
            let g = self.hol(&inp.plaq[p]);
            a += inp.ta[g];
            b += inp.tb[g];
        }
        self.psa[j] = a;
        self.psb[j] = b;
        for &o in &pl.obs_at[j] {
            let g = self.hol(&inp.obs[o]);
            self.otid[o] = inp.tid[g];
        }
    }

    fn record(&mut self) {
        let k = self.inp.k;
        let w = self.pl.width;
        let ia = (self.psa[k - 1] + self.pl.off) as usize;
        let ib = (self.psb[k - 1] + self.pl.off) as usize;
        let cell = ia * w + ib;
        self.act[cell] += 1;
        let ntr = self.inp.ntr;
        for (o, &t) in self.otid.iter().enumerate() {
            self.obsh[(o * w * w + cell) * ntr + t] += 1;
        }
    }

    fn run_task(&mut self, task: usize) {
        let n = self.inp.n;
        let k = self.inp.k;
        let pre = self.pl.prefix;
        let mut t = task;
        for i in (0..pre).rev() {
            self.digit[i] = t % n;
            t /= n;
        }
        for i in pre..k {
            self.digit[i] = 0;
        }
        for j in 0..k {
            self.eval(j);
        }
        loop {
            self.record();
            let mut i = k - 1;
            loop {
                if i < pre {
                    return;
                }
                self.digit[i] += 1;
                if self.digit[i] < n {
                    break;
                }
                self.digit[i] = 0;
                if i == 0 {
                    return;
                }
                i -= 1;
            }
            for j in i..k {
                self.eval(j);
            }
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert!(args.len() >= 3, "usage: transfer_enum <input> <output> [threads]");
    let threads: usize = if args.len() > 3 {
        args[3].parse().expect("threads")
    } else {
        8
    };
    let text = fs::read_to_string(&args[1]).expect("read input");
    let inp = parse(&text);
    assert!(inp.k >= 1, "need at least one free link");
    let pl = plan(&inp);
    let ntasks = inp.n.pow(pl.prefix as u32);
    let counter = AtomicUsize::new(0);
    let parts: Vec<(Vec<u64>, Vec<u64>)> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut wk = Walker::new(&inp, &pl);
                    loop {
                        let task = counter.fetch_add(1, Ordering::Relaxed);
                        if task >= ntasks {
                            break;
                        }
                        wk.run_task(task);
                    }
                    (wk.act, wk.obsh)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("worker"))
            .collect()
    });
    let w = pl.width;
    let mut act = vec![0u64; w * w];
    let mut obsh = vec![0u64; inp.obs.len() * w * w * inp.ntr];
    for (a, o) in &parts {
        for (x, y) in act.iter_mut().zip(a.iter()) {
            *x += *y;
        }
        for (x, y) in obsh.iter_mut().zip(o.iter()) {
            *x += *y;
        }
    }
    let total: u64 = act.iter().sum();
    let mut out = String::new();
    out.push_str(&format!("TOTAL {}\n", total));
    for ia in 0..w {
        for ib in 0..w {
            let c = act[ia * w + ib];
            if c > 0 {
                out.push_str(&format!(
                    "ACT {} {} {}\n",
                    ia as i64 - pl.off,
                    ib as i64 - pl.off,
                    c
                ));
            }
        }
    }
    for o in 0..inp.obs.len() {
        for ia in 0..w {
            for ib in 0..w {
                for t in 0..inp.ntr {
                    let c = obsh[(o * w * w + ia * w + ib) * inp.ntr + t];
                    if c > 0 {
                        out.push_str(&format!(
                            "OBS {} {} {} {} {}\n",
                            o,
                            ia as i64 - pl.off,
                            ib as i64 - pl.off,
                            t,
                            c
                        ));
                    }
                }
            }
        }
    }
    fs::write(&args[2], out).expect("write output");
}
